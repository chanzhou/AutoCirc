# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


### [1.3] - 2018-8-21
### Added
- Incoporate the function of removing false positive circRNAs that come across multiple unoverlapped transcripts, but do not originate from another transcripts;
- Include "circ_across_mgenes.pl" in the script folder.

### Changed
- Update the READ file regarding the usage of "--min" and "--max" arguments.

### Removed
-Remove "circ_across2genes.pl" and "circ2mboundaries.pl" from the script folder and related lines in the AutoCirc.pl

### [1.2] - 2017-10-12
### Fixed
- Improvement of using the standarded BED format with 12 fields for the reference gene annotation file

### Removed
- The "--sensitive" argument 


### [1.1] - 2017-06-30
### Added
- "-s/--seed" and "--sensitive" options to allow selecting different sizes of seed for initial mapping from unmapped reads.


### [1.0] - 2017-06-20
### Added
- Release v1.0
