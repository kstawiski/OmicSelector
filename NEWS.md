# OmicSelector 2.1.0 (2025-08-14)

## Major Changes

### Package Modernization
* **BREAKING CHANGE**: Minimum R version increased to 4.1.0
* Updated all dependencies to current versions
* Added proper error handling and input validation throughout
* Implemented modern R package structure with comprehensive testing

### New Features
* Added comprehensive unit testing framework with testthat
* Implemented modern logging system using `logger` package
* Added configuration management system
* Enhanced progress reporting with `cli` package
* Better integration with tidyverse workflows

### Code Quality Improvements
* Refactored main `OmicSelector_OmicSelector()` function into modular components
* Removed excessive use of `suppressMessages()` 
* Improved error messages and warnings
* Added proper input validation for all functions
* Enhanced documentation with better examples

### Machine Learning Enhancements
* Better integration with tidymodels framework
* Improved keras/tensorflow integration as optional dependencies
* Enhanced parallel processing capabilities
* Modern hyperparameter tuning approaches

### Bug Fixes
* Fixed duplicate RColorBrewer dependency in DESCRIPTION
* Resolved naming convention inconsistencies
* Fixed hardcoded path dependencies
* Improved memory management in large-scale analyses

### Documentation
* Updated all roxygen2 documentation
* Added comprehensive vignettes
* Better parameter descriptions and examples
* Improved README with modern installation instructions

---

# OmicSelector 2.0.0 (2025-08-14)

## Major Changes

* **Package Modernization**: Complete refactoring and modernization of the package
* **Updated Dependencies**: Migrated to modern R packages and removed deprecated dependencies
* **Improved Architecture**: Refactored large functions into modular, maintainable components
* **Enhanced Testing**: Added comprehensive test suite with testthat framework
* **Better Documentation**: Updated all documentation with roxygen2 7.3.2

## New Features

* **Modern Parallel Processing**: Replaced snow/doParallel with future/furrr for better performance
* **Tidymodels Integration**: Added support for modern tidymodels ecosystem
* **Enhanced Visualization**: Improved plotting capabilities with modern ggplot2 features
* **S3 Classes**: Proper object-oriented design with S3 classes for results
* **Input Validation**: Comprehensive input validation and error handling

## Feature Selection Methods

* **Updated Methods**: All 60+ feature selection methods have been updated and modernized
* **New Methods**: Added support for latest feature selection algorithms
* **Performance Improvements**: Optimized algorithms for better speed and memory usage

## Deep Learning

* **Keras Updates**: Updated Keras integration for compatibility with TensorFlow 2.x
* **GPU Support**: Enhanced GPU acceleration capabilities
* **Model Architecture**: Improved neural network architectures

## Breaking Changes

* **Function Names**: Some internal function names have been updated for consistency
* **Return Objects**: Function return values now use S3 classes instead of raw lists
* **Dependencies**: Removed deprecated packages - see migration guide in vignettes

## Bug Fixes

* Fixed memory leaks in parallel processing
* Improved error handling across all functions
* Fixed compatibility issues with latest R versions
* Resolved namespace conflicts

## Documentation

* Updated all vignettes with modern examples
* Added pkgdown website support
* Improved function documentation
* Added migration guide from v1.x

---

# OmicSelector 1.0.0 (Previous Release)

## Features

* Initial release with 60+ feature selection methods
* Benchmarking capabilities with multiple ML algorithms
* Deep learning integration with Keras
* Docker containerization
* Web-based GUI interface
* Support for miRNA-seq, RNA-seq, and proteomics data
