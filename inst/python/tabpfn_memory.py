"""Deterministic query limits for float32 cached math attention.

The budget covers three dense attention matrices, not the retained context,
parameters, or the process. It is identical on CPU and CUDA.
"""

ATTENTION_WORKSPACE_BYTES = 4 * 1024 ** 3


def query_limit(layouts, budget_bytes=ATTENTION_WORKSPACE_BYTES):
    """Return rows fitting the largest (groups, context, heads, bytes) layout."""
    budget_bytes = int(budget_bytes)
    if budget_bytes < 1:
        raise ValueError("TabPFN attention workspace budget must be positive")
    per_row = 0
    for layout in layouts:
        if len(layout) != 4 or any(int(n) != n or n < 1 for n in layout):
            raise ValueError("TabPFN attention layout must contain four positive integers")
        groups, context, heads, item_bytes = map(int, layout)
        per_row = max(per_row, 3 * groups * context * heads * item_bytes)
    if per_row > budget_bytes:
        raise MemoryError(
            "TabPFN frozen attention needs at least %d bytes for one query; "
            "the workspace budget is %d bytes. Reduce hp$max_context or the "
            "training feature count." % (per_row, budget_bytes))
    return max(1, budget_bytes // per_row) if per_row else 1024
