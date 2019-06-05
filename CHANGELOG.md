# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2019-06-04
### Added
- Change log

### Changed
- Increased the value of the symbolic INFINITY used in GMP to 1E1000

### Fixed
- Previous versions of the submodule fastpathz had an error in the scientific notation conversion whereby the weight of an edge was put into a character array that was not proprely null \0 terminated. This caused unpredictable behavior sometimes causing the weight to default to zero. 
