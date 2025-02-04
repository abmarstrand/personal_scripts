#!/usr/bin/env python3
"""
Standalone Loom File Generator for RNA Velocity

This script implements a simplified version of velocyto’s functionality:
it reads a GTF file, parses gene models (merging exon intervals),
builds interval trees for fast lookup, processes a BAM file (using pysam)
to classify reads as spliced, unspliced or ambiguous per gene, and writes a loom file
with three layers: spliced, unspliced and ambiguous.

Installation instructions:
--------------------------
This script requires Python 3.10+ and the following packages:
   - numpy
   - pysam
   - loompy
   - intervaltree

You can install the dependencies via pip:

    pip install numpy pysam loompy intervaltree

Usage:
   Adjust the paths in the __main__ section or call the function
   create_loom_from_bam_gtf(bam_file, gtf_file, output_loom, ...) with your parameters.

Note:
   • The BAM file must be coordinate-sorted by mapping position. If not, this script
     will invoke samtools sort (ensure a recent version of samtools, ≥1.6, is installed).
   • The read-classification algorithm is simplified and mimics velocyto’s default behavior:
     It assigns a read as “spliced” if all its aligned blocks fall within annotated exons,
     “unspliced” if none do (but its midpoint lies within the gene boundaries),
     and “ambiguous” otherwise.
   • If a read overlaps more than one gene, it is counted for each overlapping gene.
"""

import os
import shlex
import subprocess
import multiprocessing
import itertools
import random
import string
import logging
import gzip
import csv

import numpy as np
import pysam
import loompy
from intervaltree import Interval, IntervalTree

# -----------------------------
# Helper functions
# -----------------------------

def id_generator(size: int = 6, chars: str = string.ascii_uppercase + string.digits) -> str:
    """Generate a random ID of given size."""
    return ''.join(random.choice(chars) for _ in range(size))


def parse_gtf(gtf_file: str) -> dict:
    """
    Parse a GTF file and return a dictionary of gene models.
    
    Returns a dict mapping gene_id to a dict with keys:
      'gene_id', 'gene_name', 'chrom', 'strand', 'start', 'end', 'exons'
    where 'exons' is a sorted list of non-overlapping (merged) (start, end) tuples.
    """
    genes = {}
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, source, feature, start, end, score, strand, frame, attr = fields
            start, end = int(start), int(end)
            # Parse attributes
            attr_dict = {}
            for item in attr.strip().split(";"):
                item = item.strip()
                if item == "":
                    continue
                if " " in item:
                    key, value = item.split(" ", 1)
                    value = value.replace('"', '').strip()
                    attr_dict[key] = value
            gene_id = attr_dict.get("gene_id")
            gene_name = attr_dict.get("gene_name", gene_id)
            if gene_id is None:
                continue
            if gene_id not in genes:
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [(start, end)]
                }
            else:
                genes[gene_id]["exons"].append((start, end))
    # Merge exon intervals per gene and compute gene boundaries.
    for gene in genes.values():
        # Sort exons by start
        exons = sorted(gene["exons"], key=lambda x: x[0])
        merged = []
        current_start, current_end = exons[0]
        for s, e in exons[1:]:
            if s <= current_end:
                current_end = max(current_end, e)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = s, e
        merged.append((current_start, current_end))
        gene["exons"] = merged
        gene["start"] = merged[0][0]
        gene["end"] = merged[-1][1]
    return genes


def build_interval_trees(genes: dict) -> dict:
    """
    Build an interval tree per chromosome mapping gene regions to gene IDs.
    
    Returns a dict mapping chromosome to an IntervalTree.
    """
    trees = {}
    for gene in genes.values():
        chrom = gene["chrom"]
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        # Use gene boundaries (start, end) as the interval.
        trees[chrom][gene["start"]: gene["end"] + 1] = gene["gene_id"]
    return trees


def classify_read(read, exon_intervals: list) -> str:
    """
    Given a pysam.AlignedSegment (read) and a list of exon intervals [(s,e),...],
    classify the read as:
      - "spliced" if every aligned block lies entirely within an exon interval,
      - "unspliced" if none of the blocks overlap any exon interval,
      - "ambiguous" otherwise.
      
    Uses read.get_blocks() to retrieve aligned blocks.
    """
    blocks = read.get_blocks()  # list of (start, end) tuples
    if not blocks:
        return "ambiguous"
    contained = 0
    not_contained = 0
    for block in blocks:
        bstart, bend = block
        in_exon = any(s <= bstart and bend <= e for s, e in exon_intervals)
        if in_exon:
            contained += 1
        else:
            not_contained += 1
    if contained == len(blocks):
        return "spliced"
    elif not_contained == len(blocks):
        return "unspliced"
    else:
        return "ambiguous"


# -----------------------------
# Main processing function
# -----------------------------

def create_loom_from_bam_gtf(
    bam_file: str,
    gtf_file: str,
    output_loom: str,
    *,
    temp_dir: str = None,
    samtools_threads: int = 16,
    samtools_memory: int = 2048,
    loom_numeric_dtype: str = "uint32",
    verbose: int = 2
) -> None:
    """
    Stand-alone function to generate a loom file with spliced, unspliced and ambiguous counts.
    
    Parameters:
      bam_file (str): Path to a BAM file (may be unsorted by cell barcode).
      gtf_file (str): Path to a GTF file containing gene annotation.
      output_loom (str): Path to the output loom file.
      temp_dir (str, optional): Directory for temporary files (samtools -T flag).
      samtools_threads (int): Number of threads for samtools sort.
      samtools_memory (int): Memory (in MB) per thread for samtools sort.
      loom_numeric_dtype (str): Dtype for loom layers (e.g., "uint32").
      verbose (int): Verbosity level (0=ERROR, 1=WARNING, 2=INFO, 3=DEBUG).
      
    Returns:
      None. Writes the loom file at output_loom.
      
    The function:
      1. Parses the GTF file to build gene models.
      2. Builds per-chromosome interval trees for gene lookup.
      3. Sorts the BAM file by cell barcode (using samtools sort) if needed.
      4. Iterates through the BAM file (using pysam) to count molecules per gene per cell.
      5. For each read, determines whether it is spliced, unspliced or ambiguous.
      6. Writes the count matrices to a loom file with row attributes (gene info) and column attributes (cell IDs).
    """
    # Configure logging.
    log_levels = [logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=log_levels[verbose]
    )
    
    logging.info("Parsing GTF file...")
    genes = parse_gtf(gtf_file)
    logging.info(f"Parsed {len(genes)} genes from {gtf_file}.")
    
    logging.info("Building interval trees for gene lookup...")
    trees = build_interval_trees(genes)
    
    # Determine sorted BAM filename.
    if already_sorted:
        bam_sorted = bam_file
        logging.info("Using provided BAM file as already sorted.")
    else:
        bam_sorted = os.path.join(os.path.dirname(bam_file), f"cellsorted_{os.path.basename(bam_file)}")
        if not os.path.exists(bam_sorted):
            try:
                mem_line = subprocess.check_output(['grep', 'MemAvailable', '/proc/meminfo'])
                mb_available = int(mem_line.split()[1]) / 1000
            except Exception:
                logging.warning("Could not determine available memory; assuming 32000 MB")
                mb_available = 32000
            threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
            mb_to_use = int(min(samtools_memory, mb_available / threads_to_use))
            cmd = [
                "samtools", "sort",
                "-m", f"{mb_to_use}M",
                "-O", "BAM",
                "-@", str(threads_to_use)
            ]
            if temp_dir is not None:
                cmd += ["-T", temp_dir]
            cmd += ["-o", bam_sorted, bam_file]
            logging.info("Sorting BAM file with samtools...")
            logging.debug(f"Samtools command: {' '.join(shlex.quote(c) for c in cmd)}")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                raise MemoryError(f"Samtools sort failed (return code {proc.returncode}):\n{stderr.decode()}")
            logging.info("BAM file sorted successfully.")
        else:
            logging.info(f"Sorted BAM file exists: {bam_sorted}")
    
    # Open the sorted BAM file.
    bam_in = pysam.AlignmentFile(bam_sorted, "rb")
    
    # We assume that each read has a cell barcode in tag "CB".
    # If not present, assign a default "unknown" barcode.
    cell_counts = {}  # keys: cell barcode, values: dict with gene counts per layer
    # Structure: counts[gene_id][cell_barcode] for each layer; we create separate dictionaries.
    spliced_counts = {}
    unspliced_counts = {}
    ambiguous_counts = {}
    
    # Process reads.
    logging.info("Processing BAM file and counting reads...")
    n_processed = 0
    for read in bam_in.fetch(until_eof=True):
        # Skip unmapped or secondary/supplementary reads.
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        try:
            cell_barcode = read.get_tag("CB")
        except KeyError:
            cell_barcode = "unknown"
        # Get read midpoint.
        mid = (read.reference_start + read.reference_end) // 2
        chrom = bam_in.get_reference_name(read.reference_id)
        # Check if chromosome exists in our trees.
        if chrom not in trees:
            continue
        # Get candidate genes overlapping the midpoint.
        candidates = trees[chrom].at(mid)
        if not candidates:
            continue
        # For each candidate gene, classify the read.
        for iv in candidates:
            gene_id = iv.data
            gene_model = genes[gene_id]
            # For this gene, get its exon intervals.
            exon_intervals = gene_model["exons"]
            # Classify the read.
            classification = classify_read(read, exon_intervals)
            # Initialize nested dictionaries if needed.
            if gene_id not in spliced_counts:
                spliced_counts[gene_id] = {}
                unspliced_counts[gene_id] = {}
                ambiguous_counts[gene_id] = {}
            # Increment count for the cell.
            if classification == "spliced":
                spliced_counts[gene_id][cell_barcode] = spliced_counts[gene_id].get(cell_barcode, 0) + 1
            elif classification == "unspliced":
                unspliced_counts[gene_id][cell_barcode] = unspliced_counts[gene_id].get(cell_barcode, 0) + 1
            else:
                ambiguous_counts[gene_id][cell_barcode] = ambiguous_counts[gene_id].get(cell_barcode, 0) + 1
        n_processed += 1
        if n_processed % 100000 == 0:
            logging.info(f"Processed {n_processed} reads...")
    bam_in.close()
    logging.info(f"Finished processing BAM file. Total processed reads: {n_processed}")
    
    # Determine the set of all genes and cells observed.
    all_gene_ids = sorted(genes.keys())
    all_cell_barcodes = set()
    for d in (spliced_counts, unspliced_counts, ambiguous_counts):
        for gene in d:
            all_cell_barcodes.update(d[gene].keys())
    all_cell_barcodes = sorted(all_cell_barcodes)
    
    n_genes = len(all_gene_ids)
    n_cells = len(all_cell_barcodes)
    logging.info(f"Building count matrices: {n_genes} genes x {n_cells} cells")
    
    # Create count matrices (rows: genes, columns: cells)
    def build_matrix(counts_dict):
        mat = np.zeros((n_genes, n_cells), dtype=np.uint32)
        cell_to_ix = {cell: i for i, cell in enumerate(all_cell_barcodes)}
        gene_to_ix = {gene: i for i, gene in enumerate(all_gene_ids)}
        for gene, cell_dict in counts_dict.items():
            i = gene_to_ix.get(gene)
            if i is None:
                continue
            for cell, count in cell_dict.items():
                j = cell_to_ix[cell]
                mat[i, j] = count
        return mat
    
    spliced_matrix = build_matrix(spliced_counts)
    unspliced_matrix = build_matrix(unspliced_counts)
    ambiguous_matrix = build_matrix(ambiguous_counts)
    
    # Build row attributes (gene metadata)
    ra = {}
    ra["Gene"] = np.array([genes[gid]["gene_name"] for gid in all_gene_ids], dtype=str)
    ra["Accession"] = np.array(all_gene_ids, dtype=str)
    ra["Chromosome"] = np.array([genes[gid]["chrom"] for gid in all_gene_ids], dtype=str)
    ra["Strand"] = np.array([genes[gid]["strand"] for gid in all_gene_ids], dtype=str)
    ra["Start"] = np.array([genes[gid]["start"] for gid in all_gene_ids], dtype=int)
    ra["End"] = np.array([genes[gid]["end"] for gid in all_gene_ids], dtype=int)
    
    # Build column attributes (cell IDs)
    ca = {"CellID": np.array(all_cell_barcodes, dtype=str)}
    
    # Prepare layers dictionary
    layers = {
        "spliced": spliced_matrix.astype(loom_numeric_dtype, copy=False),
        "unspliced": unspliced_matrix.astype(loom_numeric_dtype, copy=False),
        "ambiguous": ambiguous_matrix.astype(loom_numeric_dtype, copy=False)
    }
    
    # Write loom file using loompy (compatible with loompy v2)
    logging.info("Writing loom file...")
    loompy.create(
        filename=output_loom,
        layers=layers,
        row_attrs=ra,
        col_attrs=ca,
        file_attrs={
            "script_version": "stand-alone-velocyto.py",
            "notes": "Counts were generated by a stand-alone implementation that mimics velocyto Default logic"
        }
    )
    logging.info(f"Loom file written successfully to {output_loom}.")


# -----------------------------
# Main entry point
# -----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Stand-alone loom file generator for spliced, unspliced, and ambiguous counts."
    )
    parser.add_argument("bam_file", help="Path to the BAM file.")
    parser.add_argument("gtf_file", help="Path to the GTF annotation file.")
    parser.add_argument("output_loom", help="Path to the output loom file.")
    parser.add_argument("--already_sorted", action="store_true",
                        help="Flag to indicate that the provided BAM file is already sorted.")
    parser.add_argument("--temp_dir", default=None, help="Temporary directory for samtools sort (-T flag).")
    parser.add_argument("--samtools_threads", type=int, default=16, help="Number of threads for samtools sort.")
    parser.add_argument("--samtools_memory", type=int, default=2048, help="Memory (MB) per thread for samtools sort.")
    parser.add_argument("--loom_dtype", default="uint32", help="Numeric dtype for loom layers.")
    parser.add_argument("--verbose", type=int, default=2, help="Verbosity level (0=ERROR, 1=WARNING, 2=INFO, 3=DEBUG).")
    args = parser.parse_args()
    
    create_loom_from_bam_gtf(
        bam_file=args.bam_file,
        gtf_file=args.gtf_file,
        output_loom=args.output_loom,
        already_sorted=args.already_sorted,
        temp_dir=args.temp_dir,
        samtools_threads=args.samtools_threads,
        samtools_memory=args.samtools_memory,
        loom_numeric_dtype=args.loom_dtype,
        verbose=args.verbose
    )
