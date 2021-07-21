# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.5.0] - 2020-07-09
### Fixed
- An error in the fastpath package that caused the best path not to be found if edges had large enough weights. This is because of rounding errors for long doubles, and that C cannot distinguish between 1E50 and 1E50 + 1
### Changed
- Switched from using fastpath to fastpathz, which uses infinite-precision integers

## [1.4.0] - 2020-06-10
### Fixed
- An error in the fastpath package caused by not emptying out the graph between contigs, if the input file has more than one contig

## [1.3.1] - 2020-05-25
### Added
- Missing MANIFEST.in file

## [1.3.0] - 2020-05-05
### Changed
- All changes in this version are merely optimization and organization, output should be identical to previous versions.
- Reorganized phanotate into a package that installs via setup.py script or through pip
- Coded up fastpath into a Cython pypi extension that is fetched and installed via pip
- Renamed the modules folder

## [1.2.2] - 2019-09-08
### Fixed
- An error that was caused by newer versions of tRNAscanSE using the -b flag to indicate a bed file, instead of the earlier versions where this flag was to output brief format.  Using the long arg --brief is compatible with all versions.  
- An error where the ouput from tRNAscanSE with popen was not decoded from bytes to string properly.

## [1.2.1] - 2019-06-06
### Fixed
- An error that was caused by connecting the wrong types of CDS and tRNA in the graph network

### Added
- Displaying the weights (scores) of ORFs in the genbank out file format

## [1.1.1] - 2019-06-06
### Changed
- Excluded tRNA gene calls from tabular and fasta output

### Fixed
- The reporting of tRNAs when the output format is genbank

## [1.1.0] - 2019-06-04
### Added
- Change log

### Changed
- Increased the value of the symbolic INFINITY used in GMP to 1E1000

### Fixed
- Previous versions of the submodule fastpathz had an error in the scientific notation conversion whereby the weight of an edge was put into a character array that was not proprely null \0 terminated. This caused unpredictable behavior sometimes causing the weight to default to zero. 
