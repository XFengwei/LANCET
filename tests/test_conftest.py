"""Tests for the local-input skip rule in conftest.py."""

from pathlib import Path

import pytest

from conftest import missing_local_trees


def test_missing_local_trees_names_only_absent_trees(tmp_path: Path) -> None:
    (tmp_path / "data").mkdir()
    trees = {"data": tmp_path / "data", "results": tmp_path / "results"}

    assert missing_local_trees(("data", "results"), trees) == ["results"]
    assert missing_local_trees(("data",), trees) == []


def test_unknown_tree_name_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unknown"):
        missing_local_trees(("papers",), {"data": tmp_path})
