# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-10-07
### Changed
- Replaced FASTQ-based workflow with new version that use aligned BAM + BAI files as input.
- Extracts chr6 region and alt contigs from BAM, optionally includes unmapped reads for RNA libraries.
- Reads are converted to FASTQ and filtered against HLA reference using razers3.
- Optimized razers3 execution with multithreading.
- Workflow is now more compatible with BAM-based pipelines and eliminates preprocessing steps (e.g. concatenation, slicing, read counting).
- Simplified scatter and concatenation logic by removing fastq chunking and merging.

## [1.1.0] - 2025-05-09
### Changed
- Modified the input data to allow multiple fastq per read (array).

## [Unreleased] - 2025-03-06
### Changed
- Modified WDL to make numReads optional input, similar to the bwaMem workflow. This will then run an linecount on the fastq file 

## [1.0.0] - 2024-12-10
### Added
- [GRD-831](https://jira.oicr.on.ca/browse/GRD-831), based on repo initialed last month with wdl file
- finalize wdl
- add readme and vidarr files



