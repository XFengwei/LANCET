"""Tests for project discovery and the data access policy."""

from pathlib import Path

import pytest

from lancet import DataOutputError, ProjectPaths


def test_project_root_is_discovered_from_nested_directory() -> None:
    project = ProjectPaths.discover(Path(__file__).parent)

    assert project.root == Path(__file__).resolve().parents[1]
    assert project.delivered == project.root / "data" / "delivered"
    assert project.reduction == project.root / "data" / "reduction"


def test_delivered_input_is_resolved_without_special_permission() -> None:
    project = ProjectPaths.discover()

    assert project.delivered_input("alma", "7M") == (
        project.root / "data" / "delivered" / "alma" / "7M"
    ).resolve()


def test_reduction_input_requires_explicit_declaration() -> None:
    project = ProjectPaths.discover()

    with pytest.raises(PermissionError, match="declared=True"):
        project.reduction_input("alma", "7MTM2TM1")

    assert project.reduction_input(
        "alma", "7MTM2TM1", declared=True
    ) == (project.root / "data" / "reduction" / "alma" / "7MTM2TM1").resolve()


def test_paths_cannot_escape_their_configured_root() -> None:
    project = ProjectPaths.discover()

    with pytest.raises(ValueError, match="escapes"):
        project.delivered_input("..", "reduction")


def test_outputs_inside_data_are_rejected() -> None:
    project = ProjectPaths.discover()

    with pytest.raises(DataOutputError, match="must not"):
        project.validate_output("data/delivered/generated.fits")

    assert project.validate_output("results/table.csv") == (
        project.root / "results" / "table.csv"
    )


def test_result_path_has_no_target_layer() -> None:
    project = ProjectPaths.discover()

    assert project.result("tables", "catalog.csv") == (
        project.root / "results" / "tables" / "catalog.csv"
    )


def test_guard_outputs_refuses_existing_files_without_overwrite(tmp_path: Path) -> None:
    from lancet.paths import guard_outputs

    existing = tmp_path / "old.fits"
    existing.write_bytes(b"")

    guard_outputs([tmp_path / "new.fits"], overwrite=False)
    guard_outputs([existing], overwrite=True)
    with pytest.raises(FileExistsError, match="--overwrite"):
        guard_outputs([existing, tmp_path / "new.fits"], overwrite=False)


def _project_with_linked_data(tmp_path: Path) -> ProjectPaths:
    """A checkout whose data/delivered is a symlink, as in a fresh clone."""
    root = tmp_path / "clone"
    (root / "data").mkdir(parents=True)
    (root / "pyproject.toml").write_text('[project]\nname = "lancet"\n')
    target = tmp_path / "shared" / "delivered"
    (target / "alma").mkdir(parents=True)
    (root / "data" / "delivered").symlink_to(target)
    return ProjectPaths(root)


def test_linked_delivered_data_resolves_to_its_target(tmp_path: Path) -> None:
    project = _project_with_linked_data(tmp_path)

    assert project.delivered_input("alma") == (tmp_path / "shared" / "delivered" / "alma").resolve()


def test_outputs_into_linked_data_are_rejected_as_data_outputs(tmp_path: Path) -> None:
    project = _project_with_linked_data(tmp_path)

    with pytest.raises(DataOutputError, match="must not"):
        project.validate_output("data/delivered/generated.fits")
