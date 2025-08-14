# Deep Learning Integration for OmicSelector

## Overview

The deep learning extension has been successfully integrated into the base OmicSelector package, providing comprehensive neural network capabilities for miRNA biomarker analysis.

## Integrated Functions

### Core Deep Learning Functions

1. **`create_keras_model()`** / **`OmicSelector_keras_create_model()`**
   - Creates Keras neural network models based on hyperparameters
   - Supports multiple layers, dropout, L1 regularization, and various activation functions
   - Designed for miRNA Case/Control classification

2. **`train_deep_learning()`** / **`OmicSelector_deep_learning()`**
   - Comprehensive deep learning training with hyperparameter optimization
   - Supports autoencoder for feature extraction
   - Includes parallel processing capabilities (framework ready)
   - Currently shows integration message - full implementation in progress

3. **`predict_deep_learning()`** / **`OmicSelector_deep_learning_predict()`**
   - Makes predictions using trained deep learning models
   - Handles scaling, autoencoder preprocessing, and evaluation
   - Supports both blinded and unblinded validation
   - Provides comprehensive performance metrics (AUC, accuracy, sensitivity, specificity)

4. **`transfer_learning_neural_network()`** / **`OmicSelector_transfer_learning_neural_network()`**
   - Performs transfer learning with pre-trained models
   - Supports layer freezing and fine-tuning
   - Handles scaling parameter restoration

5. **`load_deeplearning_extension()`**
   - Loads the complete original extension for full functionality
   - Bridges between integrated and original implementation

## Integration Status

### ✅ **Completed**
- Core function framework integrated into base package
- Keras model creation fully functional
- Deep learning prediction fully functional (when keras installed)
- Transfer learning framework implemented
- Backward compatibility maintained with original function names
- Comprehensive error handling and validation
- Modern R package structure with proper documentation

### 🚧 **In Progress**
- Full training pipeline integration (complex parallel processing)
- Complete autoencoder implementation
- Hyperparameter optimization workflow

### 📋 **Usage**

#### For Full Functionality (Recommended)
```r
# Use the original extension for complete functionality
OmicSelector_load_extension("deeplearning")

# Then use original functions as before
results <- OmicSelector_deep_learning(selected_miRNAs = ".")
```

#### For Integrated Functions
```r
library(OmicSelector)

# Create Keras models
hyperparams <- expand.grid(
  layer1 = c(32, 64),
  layer2 = c(16, 32),
  layer3 = 0,
  activation_function_layer1 = "relu",
  activation_function_layer2 = "relu", 
  activation_function_layer3 = "relu",
  dropout_layer1 = 0.2,
  dropout_layer2 = 0.1,
  dropout_layer3 = 0,
  layer1_regularizer = TRUE,
  layer2_regularizer = FALSE,
  layer3_regularizer = FALSE,
  optimizer = "adam"
)

model <- create_keras_model(1, hyperparams, 100)

# Make predictions with trained models
predictions <- predict_deep_learning(
  model_path = "trained_model.zip",
  new_dataset = new_data,
  blinded = FALSE
)
```

## Dependencies

### Required for Basic Functionality
- Standard OmicSelector dependencies
- All functions include appropriate error messages when optional packages missing

### Required for Full Deep Learning
- `keras` - Neural network framework
- `tensorflow` - Backend for keras
- `reticulate` - Python interface for R
- Additional packages as specified in function documentation

## Benefits of Integration

1. **Unified Interface**: Deep learning functions now part of main package
2. **Modern Architecture**: Proper error handling, validation, and documentation
3. **Backward Compatibility**: Original function names still work
4. **Gradual Migration**: Users can transition from extension to integrated functions
5. **Enhanced Documentation**: Full roxygen2 documentation for all functions
6. **Better Error Messages**: Clear guidance when dependencies missing

## Technical Details

### File Structure
- `R/deeplearning.R` - Main integration file with all deep learning functions
- `extensions/deeplearning.R` - Original extension (preserved for compatibility)
- Updated `NAMESPACE` and documentation files

### Key Improvements
- Removed 18+ `suppressMessages(library())` calls
- Modern namespace handling with explicit imports
- Comprehensive input validation
- Structured error handling with meaningful messages
- Modern logging integration
- Security improvements (removed unsafe `source()` calls)

## Future Development

The integration provides a solid foundation for future enhancements:

1. **Complete Training Pipeline**: Full integration of parallel hyperparameter optimization
2. **Enhanced Autoencoder Support**: Complete dimensionality reduction workflow
3. **Advanced Transfer Learning**: Extended pre-training and fine-tuning capabilities
4. **Model Management**: Improved model storage and versioning
5. **Visualization**: Enhanced plotting and monitoring capabilities

## Migration Guide

### For Current Users
- Continue using `OmicSelector_load_extension("deeplearning")` for full functionality
- Gradually transition to integrated functions as they become feature-complete
- Test new integrated functions alongside existing workflows

### For New Users
- Start with integrated functions for basic deep learning needs
- Use extension for advanced hyperparameter optimization and training
- Follow function documentation for dependency installation

This integration represents a significant step forward in modernizing OmicSelector's deep learning capabilities while maintaining full backward compatibility and providing a clear path for future development.
