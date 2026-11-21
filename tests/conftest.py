"""Skip tests whose local, gitignored inputs are absent (for example in a fresh clone).

Mark a test with ``@pytest.mark.local_inputs("data")`` when it reads
``data/delivered`` (or ``data/extraction``), and with
``@pytest.mark.local_inputs("data", "results")`` when it also reads products of
other tasks under ``results/``. Unmarked tests must run from the code alone.
"""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOCAL_TREES = {
    "data": PROJECT_ROOT / "data" / "delivered",
    "results": PROJECT_ROOT / "results",
}


def missing_local_trees(names, trees=LOCAL_TREES) -> list[str]:
    """Names among ``names`` whose local tree does not exist."""
    unknown = [name for name in names if name not in trees]
    if unknown:
        raise ValueError(f"local_inputs: unknown tree(s) {unknown}; use {sorted(trees)}")
    return [name for name in names if not trees[name].is_dir()]


def pytest_runtest_setup(item: pytest.Item) -> None:
    marker = item.get_closest_marker("local_inputs")
    if marker is None:
        return
    missing = missing_local_trees(marker.args)
    if missing:
        pytest.skip("local inputs not present: " + ", ".join(f"{name}/" for name in missing))
