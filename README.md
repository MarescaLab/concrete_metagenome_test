# Concrete metagenome test

This repository contains processing and analysis of a MiSeq test run of shotgun metagenomics from a concrete sample.

## Processing

1.  Read quality analysis with fastqc
2.  Adapter & quality trimming with TrimGalore!
3.  Paired end joining with FLASH
4.  Taxonomic profiling with MetaPhlan3
5.  Taxonomic profiling with SHOGUN
6.  Taxonomic profiling with Kraken

## Repository structure

-   data

    -   raw: raw data and QC metrics
    -   QC metagenomic data

-   code

    -   preprocess.Rmd drives the initial data processing
