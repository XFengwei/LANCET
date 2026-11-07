#!/usr/bin/env python
"""Run a bounded PHANGS 12m+7m+TP C18O combination test."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
import importlib.util
from pathlib import Path
import sys
import tomllib

from astropy.io import fits
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "tasks/c18o_phangs_combination/task.toml"


def load_config(path: Path) -> dict[str, object]:
    """Load the task configuration."""
    with path.open("rb") as handle:
        return tomllib.load(handle)


def project_path(value: object) -> Path:
    """Resolve a configured path relative to the repository root."""
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def resolve_paths(config: Mapping[str, object]) -> dict[str, Path]:
    """Resolve input and isolated PHANGS-test paths."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("inputs and outputs must be TOML tables")
    test_root = project_path(outputs["root"])
    paths = {name: project_path(value) for name, value in inputs.items()}
    paths.update(
        {
            "test_root": test_root,
            "key_dir": test_root / "keys",
            "imaging_root": test_root / "imaging",
            "postprocess_root": test_root / "postprocess",
            "derived_root": test_root / "derived",
            "release_root": test_root / "release",
            "vfield_root": test_root / "vfields",
            "cleanmask_root": test_root / "cleanmasks",
        }
    )
    return paths


def validate_inputs(paths: Mapping[str, Path]) -> None:
    """Check that the two MS inputs and TP FITS cube are present."""
    for name in ("ms_12m", "ms_7m"):
        path = paths[name]
        if not path.is_dir() or not (path / "table.dat").is_file():
            raise FileNotFoundError(f"CASA measurement set not found: {path}")
    tp_cube = paths["tp_cube"]
    if not tp_cube.is_file():
        raise FileNotFoundError(f"TP FITS cube not found: {tp_cube}")


def _human_size(path: Path) -> str:
    if path.is_file():
        size = path.stat().st_size
    else:
        size = sum(
            item.stat().st_size for item in path.rglob("*") if item.is_file()
        )
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}"
        size /= 1024
    raise AssertionError("unreachable")


def inspect_inputs(paths: Mapping[str, Path], config: Mapping[str, object]) -> None:
    """Report the input and pilot-window metadata without changing data."""
    validate_inputs(paths)
    for name in ("ms_12m", "ms_7m"):
        print(f"{name}: {paths[name]} ({_human_size(paths[name])})")

    with fits.open(paths["tp_cube"], memmap=True) as hdul:
        header = hdul[0].header
        data = np.squeeze(hdul[0].data)
        if data.ndim != 3:
            raise ValueError(f"Expected a 3-D TP cube, found shape {data.shape}")
        peak_channel = int(np.nanargmax(np.nanpercentile(data, 95, axis=(1, 2))))
        velocity_ms = (
            (peak_channel + 1 - float(header["CRPIX3"]))
            * float(header["CDELT3"])
            + float(header["CRVAL3"])
        )
        print(f"tp_cube: {paths['tp_cube']} ({_human_size(paths['tp_cube'])})")
        print(
            "  "
            f"shape={data.shape}, BUNIT={header.get('BUNIT')}, "
            f"SPECSYS={header.get('SPECSYS')}, "
            f"beam={float(header['BMAJ']) * 3600:.2f} arcsec"
        )
        print(f"  95th-percentile spectrum peaks at {velocity_ms / 1000:.3f} km/s")

    target = config["target"]
    spectral = config["spectral"]
    assert isinstance(target, Mapping) and isinstance(spectral, Mapping)
    nchan = round(float(target["velocity_width_kms"]) / float(spectral["channel_kms"]))
    print(
        "Pilot cube: "
        f"center={target['systemic_velocity_kms']} km/s, "
        f"width={target['velocity_width_kms']} km/s, "
        f"channel={spectral['channel_kms']} km/s (~{nchan} channels)"
    )
    print(f"Isolated output root: {paths['test_root']}")


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _clean_recipe(config: Mapping[str, object]) -> str:
    target = config["target"]
    spectral = config["spectral"]
    imaging = config["imaging"]
    assert isinstance(target, Mapping)
    assert isinstance(spectral, Mapping)
    assert isinstance(imaging, Mapping)
    imsize = [int(value) for value in imaging["imsize"]]
    cell = float(imaging["cell_arcsec"])
    robust = float(imaging["robust"])
    pblimit = float(imaging["pblimit"])
    interpolation = str(imaging.get("interpolation", "linear"))
    nchan = max(
        1,
        round(
            float(target["velocity_width_kms"])
            / float(spectral["channel_kms"])
        ),
    )
    phasecenter = f"J2000 {target['ra']} {target['dec']}"
    # PHANGS CleanCall removes literal spaces from every recipe line before
    # evaluating it. Escape them so the resulting CASA direction retains the
    # required separators.
    psf_phasecenter = (
        f"J2000\\x20{target['psf_ra']}\\x20{target['psf_dec']}"
    )
    return f'''# Pilot C18O mosaic recipe with independently controlled image and PSF centers.
vis = ""
selectdata = True
field = ""
spw = ""
timerange = ""
uvrange = ""
antenna = ""
scan = ""
observation = ""
intent = ""
datacolumn = "data"
imagename = ""
imsize = {imsize!r}
cell = ["{cell}arcsec"]
phasecenter = "{phasecenter}"
stokes = "I"
projection = "SIN"
startmodel = ""
specmode = "cube"
reffreq = ""
# Keep this explicit. analysisUtils.ValueMapping treats every four-channel SPW
# as WVR and removes it, so PHANGS cannot infer nchan for this four-channel
# pilot even though the measurement set is valid science data.
nchan = {nchan}
start = ""
width = ""
outframe = "LSRK"
veltype = "radio"
restfreq = []
interpolation = "{interpolation}"
perchanweightdensity = True
gridder = "mosaic"
facets = 1
psfphasecenter = "{psf_phasecenter}"
chanchunks = 1
wprojplanes = 1
vptable = ""
mosweight = True
aterm = True
psterm = False
wbawp = True
conjbeams = False
cfcache = ""
usepointing = False
computepastep = 360.0
rotatepastep = 360.0
pointingoffsetsigdev = []
pblimit = {pblimit}
normtype = "flatnoise"
deconvolver = "multiscale"
scales = []
nterms = 1
smallscalebias = 0.9
restoration = True
restoringbeam = "common"
pbcor = False
outlierfile = ""
weighting = "briggs"
robust = {robust}
noise = "1.0Jy"
npixels = 0
uvtaper = []
niter = 1
gain = 0.1
threshold = "0.0mJy/beam"
nsigma = 0.0
cycleniter = 200
cyclefactor = 3.0
minpsffraction = 0.5
maxpsffraction = 0.8
interactive = False
usemask = "user"
mask = ""
pbmask = 0.0
sidelobethreshold = 3.0
noisethreshold = 5.0
lownoisethreshold = 1.5
negativethreshold = 0.0
smoothfactor = 1.0
minbeamfrac = 0.3
cutthreshold = 0.01
growiterations = 75
dogrowprune = True
minpercentchange = -1.0
verbose = False
fastnoise = True
restart = True
savemodel = "none"
calcres = True
calcpsf = True
parallel = False
'''


def write_keys(paths: Mapping[str, Path], config: Mapping[str, object]) -> Path:
    """Generate the absolute-path key set required by PHANGS."""
    validate_inputs(paths)
    for name in (
        "key_dir",
        "imaging_root",
        "postprocess_root",
        "derived_root",
        "release_root",
        "vfield_root",
        "cleanmask_root",
    ):
        paths[name].mkdir(parents=True, exist_ok=True)

    target = config["target"]
    spectral = config["spectral"]
    imaging = config["imaging"]
    assert isinstance(target, Mapping)
    assert isinstance(spectral, Mapping)
    assert isinstance(imaging, Mapping)
    target_name = str(target["name"])
    product = str(spectral["product"])
    input_root = paths["ms_12m"].parent
    if paths["ms_7m"].parent != input_root or paths["tp_cube"].parent != input_root:
        raise ValueError("12m, 7m, and TP inputs must share one parent directory")

    key_dir = paths["key_dir"]
    master = f"""# Generated by scripts/c18o_phangs_combination.py
key_dir {key_dir}/
ms_root {input_root}/
imaging_root {paths['imaging_root']}/
postprocess_root {paths['postprocess_root']}/
derived_root {paths['derived_root']}/
release_root {paths['release_root']}/
singledish_root {input_root}/
cleanmask_root {paths['cleanmask_root']}/
vfield_root {paths['vfield_root']}/

ms_key ms_file_key.txt
singledish_key singledish_key.txt
cleanmask_key cleanmask_key.txt
distance_key distance_key.txt
window_key window_key.txt
config_key config_definitions.txt
target_key target_definitions.txt
moment_key moment_key.txt
derived_key derived_key.txt
imaging_key imaging_recipes.txt
linmos_key linearmosaic_definitions.txt
override_key overrides.txt
dir_key dir_key.txt
"""
    _write(key_dir / "master_key.txt", master)
    _write(
        key_dir / "target_definitions.txt",
        "# target ra dec vsys_kms width_kms\n"
        f"{target_name}\t{target['ra']}\t{target['dec']}\t"
        f"{target['systemic_velocity_kms']}\t{target['velocity_width_kms']}\n",
    )
    scales = [float(value) for value in imaging["clean_scales_arcsec"]]
    scales_literal = "[" + ",".join(str(value) for value in scales) + "]"
    _write(
        key_dir / "config_definitions.txt",
        "# type name parameters\n"
        f"line_product\t{product}\t"
        f"{{'line_tag':'c18o21','channel_kms':{float(spectral['channel_kms'])}}}\n"
        "interf_config\t12m+7m\t"
        "{'array_tags':['12m','7m']}\n"
        "interf_config\t12m+7m\t"
        f"{{'clean_scales_arcsec':{scales_literal}}}\n"
        "feather_config\t12m+7m+tp\t{'interf_config':'12m+7m'}\n",
    )
    _write(
        key_dir / "ms_file_key.txt",
        "# target project field array_tag obsnum filename\n"
        f"{target_name}\tLANCET\tall\t12m\t1\t{paths['ms_12m'].name}\n"
        f"{target_name}\tLANCET\tall\t7m\t1\t{paths['ms_7m'].name}\n",
    )
    _write(
        key_dir / "singledish_key.txt",
        "# target product filename\n"
        f"{target_name}\t{product}\t{paths['tp_cube'].name}\n",
    )
    _write(
        key_dir / "imaging_recipes.txt",
        "# config product stage recipe\n"
        f"all\t{product}\tall\tc18o_test.clean\n",
    )
    _write(key_dir / "c18o_test.clean", _clean_recipe(config))
    for filename in (
        "cleanmask_key.txt",
        "distance_key.txt",
        "window_key.txt",
        "moment_key.txt",
        "derived_key.txt",
        "linearmosaic_definitions.txt",
        "overrides.txt",
        "dir_key.txt",
    ):
        _write(key_dir / filename, "# Intentionally empty for the combination pilot.\n")
    for root_name in ("imaging_root", "postprocess_root", "derived_root", "release_root"):
        (paths[root_name] / target_name).mkdir(parents=True, exist_ok=True)
    print(f"Wrote PHANGS key set: {key_dir}")
    return key_dir / "master_key.txt"


def import_pipeline(pipeline_root: Path | None, analysis_utils: Path | None):
    """Import phangsPipeline, optionally from a local repository checkout."""
    if analysis_utils is not None:
        sys.path.insert(0, str(analysis_utils.resolve()))
    if pipeline_root is not None:
        sys.path.insert(0, str(pipeline_root.resolve()))
    if importlib.util.find_spec("phangsPipeline") is None:
        raise RuntimeError(
            "phangsPipeline is not importable. Install it or pass "
            "--pipeline-root /path/to/phangs_imaging_scripts."
        )
    import phangsPipeline as ppl

    return ppl


def configure_handler(handler, target: str, product: str, *, feather: bool = False) -> None:
    handler.set_targets(only=[target])
    handler.set_interf_configs(only=["12m+7m"])
    handler.set_line_products(only=[product])
    if feather:
        handler.set_feather_configs(only=["12m+7m+tp"])
    handler.set_no_cont_products(True)


def validate_keys(master_key: Path, ppl, target: str, product: str):
    """Parse the generated keys and confirm the intended combination."""
    handler = ppl.KeyHandler(master_key=str(master_key))
    assert target in handler.get_targets()
    assert product in handler.get_line_products()
    assert "12m+7m" in handler.get_interf_configs()
    assert "12m+7m+tp" in handler.get_feather_configs()
    assert handler.has_singledish(target=target, product=product)
    print("PHANGS KeyHandler validation passed")
    return handler


def run_stage(
    handler,
    ppl,
    target: str,
    product: str,
    overwrite: bool,
    *,
    copy_inputs: bool = True,
    clear_pointing: bool = True,
) -> None:
    """Copy/regrid/concatenate already-continuum-subtracted 12m and 7m data."""
    uv_handler = ppl.VisHandler(key_handler=handler)
    configure_handler(uv_handler, target, product)
    if copy_inputs:
        uv_handler.loop_stage_uvdata(
            do_copy=True,
            do_contsub=False,
            do_extract_line=False,
            require_full_line_coverage=True,
            statwt_line=False,
            overwrite=overwrite,
        )

    # VisHandler currently hard-codes clear_pointing=False. Override only the
    # batch extraction call so each temporary SPW has its POINTING rows cleared
    # before concatenation, then restore PHANGS's original function.
    casa_vis_routines = importlib.import_module(
        f"{ppl.__name__}.casaVisRoutines"
    )
    original_batch_extract_line = casa_vis_routines.batch_extract_line

    def batch_extract_line_with_pointing_policy(*args, **kwargs):
        kwargs["clear_pointing"] = clear_pointing
        return original_batch_extract_line(*args, **kwargs)

    if clear_pointing:
        print(
            "Line extraction will clear POINTING rows before SPW "
            "concatenation (tclean usepointing=False)."
        )
    casa_vis_routines.batch_extract_line = batch_extract_line_with_pointing_policy
    try:
        uv_handler.loop_stage_uvdata(
            do_copy=False,
            do_contsub=False,
            do_extract_line=True,
            require_full_line_coverage=True,
            statwt_line=False,
            overwrite=overwrite,
        )
    finally:
        casa_vis_routines.batch_extract_line = original_batch_extract_line


def run_imaging(
    handler,
    ppl,
    target: str,
    product: str,
    imaging: Mapping[str, object],
    cleaning: Mapping[str, object],
    *,
    dirty_only: bool,
    reuse_multiscale: bool = False,
) -> None:
    """Run joint 12m+7m imaging, either dirty-only or the PHANGS clean recipe."""
    if reuse_multiscale:
        imaging_dir = Path(handler.get_imaging_dir_for_target(target=target))
        checkpoint_root = imaging_dir / f"{target}_12m+7m_{product}_multiscale"
        required_checkpoint_products = ("image", "residual", "model", "pb", "psf")
        missing = [
            str(checkpoint_root) + "." + suffix
            for suffix in required_checkpoint_products
            if not Path(str(checkpoint_root) + "." + suffix).is_dir()
        ]
        if missing:
            raise FileNotFoundError(
                "Cannot run --stage reclean because the retained multiscale "
                "checkpoint is incomplete:\n  " + "\n  ".join(missing)
            )
        print(f"Recleaning from retained checkpoint: {checkpoint_root}")

    im_handler = ppl.ImagingHandler(key_handler=handler)
    configure_handler(im_handler, target, product)
    dynamic_sizing = bool(imaging.get("dynamic_sizing", True))
    if dirty_only:
        # PHANGS currently turns on all imaging stages when overwrite=False,
        # so overwrite=True is required here to keep this genuinely dirty-only.
        im_handler.loop_imaging(
            do_dirty_image=True,
            dynamic_sizing=dynamic_sizing,
            overwrite=True,
        )
        return

    mask_high_snr = float(cleaning["mask_high_snr"])
    mask_low_snr = float(cleaning["mask_low_snr"])
    mask_absolute = bool(cleaning["mask_absolute"])
    mask_spectral_roll = bool(cleaning["mask_spectral_roll"])
    mask_pb_limit = float(cleaning["mask_pb_limit"])
    singlescale_threshold_snr = float(cleaning["singlescale_threshold_snr"])
    convergence_fracflux = float(cleaning["convergence_fracflux"])
    if not 0.0 < mask_low_snr <= mask_high_snr:
        raise ValueError("Require 0 < mask_low_snr <= mask_high_snr")
    if not 0.0 < mask_pb_limit < 1.0:
        raise ValueError("mask_pb_limit must lie between zero and one")
    if singlescale_threshold_snr <= 0.0:
        raise ValueError("singlescale_threshold_snr must be positive")

    # PHANGS's large-cube mask workaround replaces a degree-based FITS WCS
    # with the radian values returned by CASA imhead, but leaves CUNIT='deg'.
    # This pilot mask is only about 32 MiB, so use CASA's direct putchunk path
    # and preserve the image coordinate system exactly.
    masking_routines = importlib.import_module(
        f"{ppl.__name__}.casaMaskingRoutines"
    )
    original_write_mask = masking_routines.write_mask
    original_signal_mask = masking_routines.signal_mask

    def write_mask_without_fits_wcs_roundtrip(*args, **kwargs):
        kwargs["huge_cube_workaround"] = False
        return original_write_mask(*args, **kwargs)

    def signal_mask_with_pilot_policy(*args, **kwargs):
        # np.roll wraps the first and last planes together. That is not a
        # physically meaningful spectral neighbour relation for this bounded
        # four-channel pilot, so make the policy configurable and default it
        # off. PHANGS writes to cube_root + '.mask' + suffix_out even though it
        # also accepts out_file; derive the actual output name the same way.
        kwargs["do_roll"] = mask_spectral_roll
        result = original_signal_mask(*args, **kwargs)

        cube_root = kwargs.get("cube_root")
        suffix_in = kwargs.get("suffix_in", "")
        suffix_out = kwargs.get("suffix_out", "")
        if cube_root is None:
            raise ValueError("signal_mask requires cube_root")
        if kwargs.get("imaging_method", "tclean") == "sdintimaging":
            cube_root += ".joint.cube"
        mask_file = cube_root + ".mask" + suffix_out
        pb_file = cube_root + ".pb" + suffix_in
        if not Path(pb_file).is_dir():
            raise FileNotFoundError(
                f"Cannot apply PB-limited clean mask; PB image not found: {pb_file}"
            )
        masking_routines.apply_additional_mask(
            old_mask_file=mask_file,
            new_mask_file=pb_file,
            new_thresh=mask_pb_limit,
            operation="AND",
        )
        print(
            f"Restricted the single-scale signal mask to PB > {mask_pb_limit} "
            f"(spectral roll={mask_spectral_roll})."
        )
        return result

    masking_routines.write_mask = write_mask_without_fits_wcs_roundtrip
    masking_routines.signal_mask = signal_mask_with_pilot_policy
    try:
        im_handler.loop_imaging(
            imaging_method="tclean",
            do_dirty_image=not reuse_multiscale,
            do_revert_to_dirty=not reuse_multiscale,
            do_read_clean_mask=not reuse_multiscale,
            do_multiscale_clean=not reuse_multiscale,
            do_revert_to_multiscale=True,
            do_singlescale_mask=True,
            singlescale_mask_absolute=mask_absolute,
            singlescale_mask_high_snr=mask_high_snr,
            singlescale_mask_low_snr=mask_low_snr,
            do_singlescale_clean=True,
            do_revert_to_singlescale=True,
            convergence_fracflux=convergence_fracflux,
            singlescale_threshold_value=singlescale_threshold_snr,
            do_export_to_fits=True,
            dynamic_sizing=dynamic_sizing,
            overwrite=True,
        )
    finally:
        masking_routines.signal_mask = original_signal_mask
        masking_routines.write_mask = original_write_mask


def run_postprocess(handler, ppl, target: str, product: str) -> None:
    """Primary-beam correct, feather the TP cube, and export final products."""
    pp_handler = ppl.PostProcessHandler(key_handler=handler)
    configure_handler(pp_handler, target, product, feather=True)
    pp_handler.loop_postprocess(
        do_prep=True,
        do_feather=True,
        do_mosaic=False,
        do_cleanup=True,
        feather_noapod=True,
        imaging_method="tclean",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=(
            "inspect",
            "keys",
            "validate",
            "stage",
            "extract",
            "dirty",
            "image",
            "reclean",
            "postprocess",
            "all",
        ),
        default="inspect",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--pipeline-root", type=Path)
    parser.add_argument("--analysis-utils", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config.resolve())
    paths = resolve_paths(config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
        return

    master_key = write_keys(paths, config)
    if args.stage == "keys":
        return

    ppl = import_pipeline(args.pipeline_root, args.analysis_utils)
    target = config["target"]
    spectral = config["spectral"]
    staging = config["staging"]
    imaging = config["imaging"]
    cleaning = config["cleaning"]
    assert isinstance(target, Mapping)
    assert isinstance(spectral, Mapping)
    assert isinstance(staging, Mapping)
    assert isinstance(imaging, Mapping)
    assert isinstance(cleaning, Mapping)
    target_name = str(target["name"])
    product = str(spectral["product"])
    key_handler = validate_keys(master_key, ppl, target_name, product)
    if args.stage == "validate":
        return

    if importlib.util.find_spec("casatasks") is None:
        raise RuntimeError(
            "CASA tasks are unavailable. Run this stage in Python 3.12 with "
            "CASA >=6.7.3 and the PHANGS CASA extras installed."
        )
    if args.stage in ("stage", "extract", "all"):
        run_stage(
            key_handler,
            ppl,
            target_name,
            product,
            args.overwrite,
            copy_inputs=args.stage != "extract",
            clear_pointing=bool(staging.get("clear_pointing", True)),
        )
    if args.stage in ("dirty", "all"):
        run_imaging(
            key_handler,
            ppl,
            target_name,
            product,
            imaging,
            cleaning,
            dirty_only=args.stage == "dirty",
        )
    if args.stage in ("image", "reclean"):
        run_imaging(
            key_handler,
            ppl,
            target_name,
            product,
            imaging,
            cleaning,
            dirty_only=False,
            reuse_multiscale=args.stage == "reclean",
        )
    if args.stage in ("postprocess", "all"):
        run_postprocess(key_handler, ppl, target_name, product)


if __name__ == "__main__":
    main()
