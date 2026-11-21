"""Tests for the gas_to_dust_ratio task."""

from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "gas_to_dust_ratio_script", PROJECT_ROOT / "tasks" / "gas_to_dust_ratio" / "run.py"
)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)
CONFIG = PROJECT_ROOT / "tasks" / "gas_to_dust_ratio" / "task.toml"


def test_config_uses_the_paper_radius() -> None:
    config = TASK.load_task_config(CONFIG)

    assert config["galaxy"]["galactocentric_radius_kpc"] == 6.5
    assert config["galaxy"]["galactocentric_radius_error_kpc"] == 0.2
    assert config["outputs"]["root"] == "results/gas_to_dust_ratio"


def test_estimate_writes_the_ratio(tmp_path: Path) -> None:
    config = TASK.load_task_config(CONFIG)

    output = TASK.build_estimate(config, tmp_path)

    result = json.loads(output.read_text())
    assert result["ratio"] == pytest.approx(10**2.0055)
    assert result["radius_kpc"] == 6.5
    assert "Giannetti" in result["relation"]


def test_existing_estimate_needs_overwrite(tmp_path: Path) -> None:
    config = TASK.load_task_config(CONFIG)
    TASK.build_estimate(config, tmp_path)

    with pytest.raises(FileExistsError):
        TASK.build_estimate(config, tmp_path)
