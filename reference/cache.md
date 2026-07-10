# Split-Aware Caching for OmicSelector

Caching utilities that prevent data leakage by including split
identifiers in cache keys. Uses memoise with cachem for efficient LRU
caching.

## Details

CRITICAL for zero-leakage: Cache keys MUST include: - Split/fold
identifier (prevents cross-fold contamination) - Preprocessing
parameters (ensures correct pipeline state) - Feature subset hash (for
filter-specific caching)

This module ensures that cached filter scores or model artifacts from
one CV fold are never reused in another fold.
