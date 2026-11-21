"""Tests for the freefree_core_contamination task."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "freefree_core_contamination_script",
    PROJECT_ROOT / "tasks" / "freefree_core_contamination" / "run.py",
)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)
CONFIG = PROJECT_ROOT / "tasks" / "freefree_core_contamination" / "task.toml"


def test_spectral_factor_is_optically_thin_scaling() -> None:
    settings = {"frequency_ghz": 223.8, "radio_frequency_ghz": 5.5, "spectral_index": -0.1}

    assert TASK.spectral_factor(settings) == pytest.approx((223.8 / 5.5) ** -0.1)


def test_footprints_are_read_by_getsf_number(tmp_path: Path) -> None:
    region = tmp_path / "foot.reg"
    region.write_text(
        "# Region file format: DS9\nfk5\n"
        "# ellipse(221.3,-59.8,1\",0.5\",30) # 7\n"
        "# ellipse(221.31,-59.81,2\",1\",0) # 12\n"
    )

    footprints = TASK.read_footprints(region)

    assert sorted(footprints) == [7, 12]


def test_inputs_are_delivered_or_task_products() -> None:
    inputs = TASK.load_task_config(CONFIG)["inputs"]

    assert inputs["continuum_image"].endswith("cont.clean.truesky.fits")
    assert inputs["aca_image"].endswith("IRDC316.8.7M.robust0.5.cont.image.fits")
    for value in inputs.values():
        assert value.startswith(("data/delivered/", "data/extraction/", "results/"))
        assert "casa" not in value and "drogdeg" not in value


def test_existing_products_need_overwrite(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text("{}")
    paths = {"output_root": tmp_path}

    with pytest.raises(FileExistsError):
        TASK.build_contamination(paths, TASK.load_task_config(CONFIG))
