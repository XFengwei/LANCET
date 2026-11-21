"""Tests for the vvv_mosaic task."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

import pytest

from lancet import ProjectPaths


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "vvv_mosaic_script", PROJECT_ROOT / "tasks" / "vvv_mosaic" / "run.py"
)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)
CONFIG = PROJECT_ROOT / "tasks" / "vvv_mosaic" / "task.toml"


def test_tiles_resolve_to_four_delivered_files_per_band() -> None:
    config = TASK.load_task_config(CONFIG)
    project = ProjectPaths(PROJECT_ROOT)
    if not (project.root / config["inputs"]["tile_directory"]).is_dir():
        pytest.skip("VVV tiles are not available")

    paths = TASK.tile_paths(project, config)

    assert set(paths) == {"j", "h", "ks"}
    assert all(len(files) == 4 for files in paths.values())


def test_tile_directory_must_be_delivered(tmp_path: Path) -> None:
    config = TASK.load_task_config(CONFIG)
    config["inputs"]["tile_directory"] = "results/vvv_mosaic"

    with pytest.raises(ValueError, match="delivered"):
        TASK.tile_paths(ProjectPaths(PROJECT_ROOT), config)


def test_existing_outputs_are_listed_per_stage(tmp_path: Path) -> None:
    config = TASK.load_task_config(CONFIG)
    (tmp_path / "vvv_ks_mosaic.fits").write_bytes(b"")
    png = tmp_path / str(config["figures"]["png"])
    png.write_bytes(b"")

    assert TASK.existing_outputs(tmp_path, config, "mosaic") == [tmp_path / "vvv_ks_mosaic.fits"]
    assert TASK.existing_outputs(tmp_path, config, "rgb") == [png]
    assert TASK.parse_args(["--overwrite"]).overwrite
    assert not TASK.parse_args([]).overwrite
