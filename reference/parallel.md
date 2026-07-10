# Parallelization Support for OmicSelector

Functions for configuring and using future-based parallelization. Uses
the \`future\` package for backend-agnostic parallel execution.

## Details

OmicSelector uses \`future\` for parallelization because: - Works with
local cores, HPC clusters, and cloud backends - Integrates natively with
mlr3 ecosystem - Supports lazy evaluation and proper error handling
