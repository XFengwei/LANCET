"""Tests for the combine_continuum task: configuration, input provenance, stage selection."""

import os
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys

import pytest

from lancet import ProjectPaths

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "tasks" / "combine_continuum" / "run.py"
SPEC = spec_from_file_location("combine_continuum_run", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "combine_continuum" / "task.toml"
RESULTS = PROJECT_ROOT / "results" / "combine_continuum"


def _config(changes):
    config = TASK.load_task_config(CONFIG)
    for (table, key), value in changes.items():
        config[table][key] = value
    return config


def _resolve(config):
    return TASK.resolve_task_paths(ProjectPaths.discover(PROJECT_ROOT), config)


@pytest.mark.parametrize("name", TASK.DELIVERED_INPUTS)
def test_delivered_inputs_must_be_delivered(name):
    config = _config({("inputs", name): "results/synthesize_continuum/model/dust_224GHz.fits"})
    with pytest.raises(ValueError, match=name):
        _resolve(config)


def test_tp_inputs_must_come_from_synthesize_continuum():
    config = _config({("inputs", "tp_dust"): "data/delivered/sedfitter/Dustmodel_223.8GHz.fits"})
    with pytest.raises(ValueError, match="tp_dust"):
        _resolve(config)


def test_measurement_sets_must_be_reduction_inputs():
    config = _config({("measurement_sets", "aca_7m"): "data/delivered/alma/continuum"})
    with pytest.raises(ValueError):
        _resolve(config)


def test_config_points_at_decided_inputs():
    config = TASK.load_task_config(CONFIG)
    assert config["inputs"]["interferometer"].endswith("cont.clean.truesky.pb0.2.fits")
    assert config["inputs"]["tp_dust_ff"] == "results/synthesize_continuum/model/continuum_224GHz.fits"
    assert "tp_frequency_ghz" not in config["prepare"]
    assert "psd_matched" not in config["tp_factor"]
    assert config["combine"]["fill_outside_mosaic"] is False
    assert "tp_display_scale" not in config["figure2"]


def test_resolved_groups():
    paths = _resolve(TASK.load_task_config(CONFIG))
    assert paths["inputs"]["gain"] == paths["gain"]
    assert paths["products"]["jcomb_continuum"].name == "jcomb_continuum_224GHz.fits"
    assert paths["diagnostics_root"] == paths["output_root"] / "diagnostics"
    assert set(paths["measurement_sets"]) == {"combined_12m", "configuration_c5", "configuration_c2", "aca_7m"}
    assert paths["regions"]["successful_combinations"].name == "success.reg"


def test_product_paths_follow_label(tmp_path):
    products = TASK.product_paths(tmp_path, "250GHz")
    assert set(products) == {"jcomb_continuum", "jcomb_dust", "feather_continuum", "feather_dust"}
    assert products["feather_dust"] == tmp_path / "products" / "feather_dust_250GHz.fits"


def test_stage_selection():
    assert TASK.stages_for("combine") == list(TASK.COMBINE_STAGES)
    assert TASK.stages_for("diagnostics") == list(TASK.DIAGNOSTIC_STAGES)
    assert TASK.stages_for("all") == [*TASK.COMBINE_STAGES, *TASK.DIAGNOSTIC_STAGES]
    assert TASK.stages_for("jcomb") == ["jcomb"]
    with pytest.raises(ValueError):
        TASK.stages_for("nope")


def test_require_outputs_names_the_stage(tmp_path):
    with pytest.raises(FileNotFoundError, match="jcomb"):
        TASK.require_outputs([tmp_path / "x.fits"], "jcomb")
    (tmp_path / "x.fits").write_text("x")
    TASK.require_outputs([tmp_path / "x.fits"], "jcomb")


def _files():
    return {p for p in RESULTS.rglob("*") if p.is_file()} if RESULTS.exists() else set()


def test_inspect_writes_nothing():
    before = _files()
    run = subprocess.run(
        [sys.executable, "tasks/combine_continuum/run.py", "--stage", "inspect"],
        cwd=PROJECT_ROOT, env={**os.environ, "PYTHONPATH": "src"}, capture_output=True, text=True,
    )
    assert run.returncode == 0, run.stderr
    assert "TP frequency: 223.8 GHz" in run.stdout
    assert "truesky.pb0.2" in run.stdout
    assert _files() == before
