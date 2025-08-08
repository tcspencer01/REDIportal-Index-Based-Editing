# REDIportal-Index-Based-Editing

# Calculate the Alu Editing Index (AEI) using REDIportal Sites

The Python script in this repository calculates the **Alu Editing Index (AEI)** from RNA-seq data by leveraging known A-to-I RNA editing sites cataloged in the **REDIportal** database. The script processes read alignments from a BAM file, isolates regions of interest, and computes a strand-specific RNA editing index based on observed A-to-G mismatches.

### What is RNA Editing?
**RNA editing** is a post-transcriptional modification that alters RNA sequences after transcription. The most common type in humans is **A-to-I editing**, where adenosine (A) is deaminated to inosine (I), which is interpreted as guanosine (G) during sequencing. A-to-I editing plays roles in neuronal function, development, and disease.

### What is the Alu Editing Index (AEI)?
The **Alu Editing Index** is a measure of A-to-I editing computed across regions enriched for editing events [1]. In this analysis, an AEI value between 0 and 1 is calculated for each region of interest in each sample/individual.

### Function
This script:

1. Accepts a BAM file of RNA-seq alignments and a BED file of REDIportal sites.
2. Accepts a BED file of **regions of interest** (e.g., Alu elements or other dsRNA regions).
3. Iterates over each region, calculating:
   - The number of A and G (or T and C) bases at REDIportal sites
   - The editing index: **G / (A + G)** on the plus strand, or **C / (T + C)** on the minus strand
4. Outputs a tab-delimited file with the AEI and supporting stats for each region

---

## Input Files

- `--input_bam (-b)`: BAM file of aligned RNA-seq reads
- `--rediportal (-r)`: BED file of REDIportal A-to-I editing sites
- `--genomic_regions (-g)`: BED file of regions to analyze, with strand encoded in the region name
- `--chr (-c)`: Chromosome to analyze (e.g., "chr7")
- `--coverage_threshold (-t)`: Minimum total coverage across REDIportal sites in a region (default = 0)
- `--output_suffix (-o)`: Optional suffix for output file name (default = `.editing_levels.txt`)

---

## Output

The output is a tab-delimited file with the following columns:
- `region_name`: Identifier for the region analyzed
- `editing_index`: AEI value for the region
- `total_REDI_sites`: Number of REDIportal sites in the region
- `covered_REDI_sites`: Number of REDIportal sites with coverage
- `REDI_sites_total_coverage`: Total coverage across all REDI sites in the region
- `total_reads`: Sum of reads across the region
- `total_positions_in_region`: Region length
- `output_strand`: Strand used for AEI calculation

## References
[1] https://pubmed.ncbi.nlm.nih.gov/31636457/

---
