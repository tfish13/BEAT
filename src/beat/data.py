"""Input adapters for cubes and collections of one-dimensional spectra."""

from __future__ import annotations

import glob
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from .adapters import AdapterError, inspect_adapter_hdul
from .config import ConfigError
from .spectrum import Spectrum


_WAVELENGTH_FACTORS = {
    "angstrom": 1.0,
    "angstroms": 1.0,
    "aa": 1.0,
    "nm": 10.0,
    "micron": 1.0e4,
    "microns": 1.0e4,
    "um": 1.0e4,
    "m": 1.0e10,
}


def wavelength_factor(unit: str) -> float:
    try:
        return _WAVELENGTH_FACTORS[str(unit).lower()]
    except KeyError as exc:
        choices = ", ".join(sorted(_WAVELENGTH_FACTORS))
        raise ConfigError(f"Unknown wavelength unit {unit!r}; choose from {choices}") from exc


def linear_wavelength(
    header: Any,
    n_pixels: int,
    fits_axis: int,
    unit: str = "angstrom",
    log10: bool = False,
) -> np.ndarray:
    """Build a FITS linear wavelength axis, including the CRPIX offset."""
    crval_key = f"CRVAL{fits_axis}"
    crpix_key = f"CRPIX{fits_axis}"
    cdelt_key = f"CDELT{fits_axis}"
    cd_key = f"CD{fits_axis}_{fits_axis}"
    if crval_key not in header:
        raise ConfigError(f"Wavelength header is missing {crval_key}")
    if cdelt_key in header:
        delta = float(header[cdelt_key])
    elif cd_key in header:
        delta = float(header[cd_key])
    else:
        raise ConfigError(f"Wavelength header is missing {cdelt_key} or {cd_key}")

    crval = float(header[crval_key])
    crpix = float(header.get(crpix_key, 1.0))
    pixels_fits = np.arange(n_pixels, dtype=float) + 1.0
    values = crval + (pixels_fits - crpix) * delta
    if log10:
        values = np.power(10.0, values)
    return values * wavelength_factor(unit)


def _uncertainty(
    values: np.ndarray | None,
    kind: str,
    scale: float,
) -> np.ndarray | None:
    if values is None:
        return None
    values = np.asarray(values, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        if kind == "sigma":
            result = np.abs(values) * scale
        elif kind == "variance":
            result = np.sqrt(values) * scale
        elif kind == "inverse_variance":
            result = np.where(values > 0, 1.0 / np.sqrt(values), np.nan) * scale
        else:
            raise ConfigError(
                "uncertainty_kind must be sigma, variance, or inverse_variance"
            )
    return result


def _slice_range(value: Any, size: int, name: str) -> range:
    if value is None:
        return range(size)
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ConfigError(f"{name} must be [start, stop] with a half-open stop")
    start = 0 if value[0] is None else int(value[0])
    stop = size if value[1] is None else int(value[1])
    if start < 0 or stop > size or start >= stop:
        raise ConfigError(f"{name}={value!r} is outside 0:{size}")
    return range(start, stop)


def _fits() -> Any:
    try:
        from astropy.io import fits
    except ImportError as exc:  # pragma: no cover - installation error
        raise RuntimeError("Astropy is required to read FITS inputs") from exc
    return fits


def _wavelength_from_hdul(
    hdul: Any,
    config: dict[str, Any],
    n_pixels: int,
    default_header_hdu: int | str,
    fits_axis: int,
) -> np.ndarray:
    unit = config.get("wavelength_unit", "angstrom")
    if config.get("adapter") in {"nirspec", "miri"} and "WCS-TABLE" in hdul:
        table = hdul["WCS-TABLE"].data
        if table is None or "wavelength" not in table.names:
            raise ConfigError("NIRSpec WCS-TABLE is missing its wavelength column")
        wavelength = np.asarray(table["wavelength"][0], dtype=float).reshape(-1)
        wavelength = wavelength * wavelength_factor("micron")
    elif "wavelength_hdu" in config:
        wavelength = np.asarray(hdul[config["wavelength_hdu"]].data, dtype=float).squeeze()
        if wavelength.ndim != 1:
            raise ConfigError("wavelength_hdu must contain a one-dimensional array")
        wavelength = wavelength * wavelength_factor(unit)
    else:
        header_hdu = config.get("wavelength_header_hdu", default_header_hdu)
        wavelength = linear_wavelength(
            hdul[header_hdu].header,
            n_pixels,
            fits_axis=int(config.get("wavelength_fits_axis", fits_axis)),
            unit=unit,
            log10=bool(config.get("wavelength_log10", False)),
        )
    if wavelength.size != n_pixels:
        raise ConfigError(
            f"Wavelength axis has {wavelength.size} pixels, expected {n_pixels}"
        )
    return wavelength


def iter_cube(config: dict[str, Any]) -> Iterator[Spectrum]:
    """Yield spectra from a three-dimensional FITS cube."""
    fits = _fits()
    path = Path(config["path"])
    with fits.open(
        path,
        memmap=True,
        do_not_scale_image_data=bool(config.get("raw_integer_mask", False)),
    ) as hdul:
        try:
            adapter_metadata = inspect_adapter_hdul(hdul, config)
        except AdapterError as exc:
            raise ConfigError(str(exc)) from exc
        flux_hdu = config["flux_hdu"]
        flux = np.asarray(hdul[flux_hdu].data)
        if flux.ndim != 3:
            raise ConfigError(f"Cube flux HDU must be 3-D; got shape {flux.shape}")
        spectral_axis = int(config.get("spectral_axis", 0))
        if spectral_axis not in {0, 1, 2}:
            raise ConfigError("input.spectral_axis must be 0, 1, or 2")
        flux = np.moveaxis(flux, spectral_axis, 0)

        uncertainty_cube = None
        if config.get("uncertainty_hdu") is not None:
            uncertainty_cube = np.asarray(hdul[config["uncertainty_hdu"]].data)
            uncertainty_cube = np.moveaxis(uncertainty_cube, spectral_axis, 0)
            if uncertainty_cube.shape != flux.shape:
                raise ConfigError("Cube flux and uncertainty HDUs have different shapes")

        mask_cube = None
        spatial_mask = None
        if config.get("mask_hdu") is not None:
            raw_mask = np.asarray(hdul[config["mask_hdu"]].data)
            if raw_mask.ndim == 3:
                mask_cube = np.moveaxis(raw_mask, spectral_axis, 0)
                if mask_cube.shape != flux.shape:
                    raise ConfigError("3-D mask and flux HDUs have different shapes")
            elif raw_mask.ndim == 2:
                spatial_mask = raw_mask.astype(bool)
                if spatial_mask.shape != flux.shape[1:]:
                    raise ConfigError("2-D mask and cube spatial shapes differ")
            else:
                raise ConfigError("mask_hdu must be a 2-D spatial or 3-D pixel mask")

        n_wave, ny, nx = flux.shape
        fits_axis = 3 - spectral_axis
        wavelength = _wavelength_from_hdul(
            hdul, config, n_wave, flux_hdu, fits_axis=fits_axis
        )
        flux_scale = float(config.get("flux_scale", 1.0))
        uncertainty_scale = float(config.get("uncertainty_scale", flux_scale))
        uncertainty_kind = config.get("uncertainty_kind", "sigma")
        x_values = _slice_range(config.get("x_range"), nx, "input.x_range")
        y_values = _slice_range(config.get("y_range"), ny, "input.y_range")
        target = config.get("target_id")
        if target is None and config.get("target_header"):
            target = hdul[config.get("target_header_hdu", 0)].header.get(
                config["target_header"]
            )
        target = str(target or path.stem).strip()
        redshift = float(config["redshift"])

        for y in y_values:
            for x in x_values:
                if spatial_mask is not None and spatial_mask[y, x]:
                    continue
                raw_uncertainty = (
                    None if uncertainty_cube is None else uncertainty_cube[:, y, x]
                )
                pixel_mask = None
                if mask_cube is not None:
                    raw_pixel_mask = np.asarray(mask_cube[:, y, x])
                    if config.get("mask_bits") is None:
                        pixel_mask = raw_pixel_mask.astype(bool)
                    else:
                        pixel_mask = (
                            np.bitwise_and(raw_pixel_mask, int(config["mask_bits"])) != 0
                        )
                yield Spectrum(
                    spectrum_id=f"{target}_x{x:04d}_y{y:04d}",
                    wavelength=wavelength,
                    flux=np.asarray(flux[:, y, x], dtype=float) * flux_scale,
                    uncertainty=_uncertainty(
                        raw_uncertainty, uncertainty_kind, uncertainty_scale
                    ),
                    redshift=redshift,
                    mask=pixel_mask,
                    metadata={
                        **adapter_metadata,
                        "input_file": str(path),
                        "x": x,
                        "y": y,
                        "fits_x": x + 1,
                        "fits_y": y + 1,
                    },
                )


def _cell_array(row: Any, column: str) -> np.ndarray:
    # Detach variable-length FITS cells from the memory-mapped table before a
    # Spectrum is handed to another process.
    return np.asarray(row[column]).squeeze().reshape(-1).copy()


def _text_id(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    return str(value).strip()


def iter_survey_table(config: dict[str, Any]) -> Iterator[Spectrum]:
    """Yield one vector spectrum per row of a FITS binary table."""
    fits = _fits()
    path = Path(config["path"])
    with fits.open(path, memmap=True) as hdul:
        table = hdul[config["hdu"]].data
        if table is None:
            raise ConfigError("Survey table HDU contains no data")
        start = int(config.get("row_start", 0))
        stop = int(config.get("row_stop", len(table)))
        if start < 0 or stop > len(table) or start >= stop:
            raise ConfigError(f"Survey row range [{start}, {stop}) is invalid")
        wave_scale = wavelength_factor(config.get("wavelength_unit", "angstrom"))
        flux_scale = float(config.get("flux_scale", 1.0))
        uncertainty_scale = float(config.get("uncertainty_scale", flux_scale))
        uncertainty_kind = config.get("uncertainty_kind", "sigma")

        for row_index in range(start, stop):
            row = table[row_index]
            uncertainty = None
            if config.get("uncertainty_column"):
                uncertainty = _uncertainty(
                    _cell_array(row, config["uncertainty_column"]),
                    uncertainty_kind,
                    uncertainty_scale,
                )
            redshift = (
                float(row[config["redshift_column"]])
                if config.get("redshift_column")
                else float(config["redshift"])
            )
            mask = None
            if config.get("mask_column"):
                mask = _cell_array(row, config["mask_column"]).astype(bool)
            yield Spectrum(
                spectrum_id=_text_id(row[config["id_column"]]),
                wavelength=_cell_array(row, config["wavelength_column"]) * wave_scale,
                flux=_cell_array(row, config["flux_column"]) * flux_scale,
                uncertainty=uncertainty,
                redshift=redshift,
                mask=mask,
                metadata={"input_file": str(path), "row": row_index},
            )


def _resolve_file_list(config: dict[str, Any]) -> list[Path]:
    if config.get("files"):
        base = Path(config.get("_config_dir", "."))
        paths = [Path(item).expanduser() for item in config["files"]]
        paths = [path if path.is_absolute() else base / path for path in paths]
    else:
        paths = [Path(item) for item in glob.glob(config["glob"], recursive=True)]
    paths = sorted(path.resolve() for path in paths if path.is_file())
    if not paths:
        raise ConfigError("spectrum_files input matched no files")
    return paths


def _file_redshift(config: dict[str, Any], header: Any | None = None) -> float:
    if config.get("redshift_header"):
        if header is None or config["redshift_header"] not in header:
            raise ConfigError(f"Missing redshift header {config['redshift_header']}")
        return float(header[config["redshift_header"]])
    if "redshift" in config:
        return float(config["redshift"])
    raise ConfigError("Each spectrum file needs redshift or redshift_header")


def _iter_fits_table_file(path: Path, config: dict[str, Any]) -> Iterator[Spectrum]:
    fits = _fits()
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul[config.get("hdu", 1)]
        table = hdu.data
        wave = np.asarray(table[config["wavelength_column"]]).squeeze().reshape(-1)
        flux = np.asarray(table[config["flux_column"]]).squeeze().reshape(-1)
        raw_uncertainty = None
        if config.get("uncertainty_column"):
            raw_uncertainty = np.asarray(
                table[config["uncertainty_column"]]
            ).squeeze().reshape(-1)
        redshift = _file_redshift(config, hdu.header)
        identifier = path.stem
        if config.get("id_header"):
            identifier = _text_id(hdu.header[config["id_header"]])
        yield Spectrum(
            spectrum_id=identifier,
            wavelength=wave * wavelength_factor(config.get("wavelength_unit", "angstrom")),
            flux=flux * float(config.get("flux_scale", 1.0)),
            uncertainty=_uncertainty(
                raw_uncertainty,
                config.get("uncertainty_kind", "sigma"),
                float(config.get("uncertainty_scale", config.get("flux_scale", 1.0))),
            ),
            redshift=redshift,
            metadata={"input_file": str(path)},
        )


def _iter_fits_image_file(path: Path, config: dict[str, Any]) -> Iterator[Spectrum]:
    fits = _fits()
    with fits.open(path, memmap=True) as hdul:
        flux_hdu = config.get("flux_hdu", 0)
        flux = np.asarray(hdul[flux_hdu].data, dtype=float).squeeze()
        if flux.ndim != 1:
            raise ConfigError(f"{path}: image spectrum is not one-dimensional")
        raw_uncertainty = None
        if config.get("uncertainty_hdu") is not None:
            raw_uncertainty = np.asarray(
                hdul[config["uncertainty_hdu"]].data, dtype=float
            ).squeeze()
        wave = _wavelength_from_hdul(
            hdul, config, flux.size, flux_hdu, fits_axis=1
        )
        header = hdul[flux_hdu].header
        identifier = path.stem
        if config.get("id_header"):
            identifier = _text_id(header[config["id_header"]])
        yield Spectrum(
            spectrum_id=identifier,
            wavelength=wave,
            flux=flux * float(config.get("flux_scale", 1.0)),
            uncertainty=_uncertainty(
                raw_uncertainty,
                config.get("uncertainty_kind", "sigma"),
                float(config.get("uncertainty_scale", config.get("flux_scale", 1.0))),
            ),
            redshift=_file_redshift(config, header),
            metadata={"input_file": str(path)},
        )


def _iter_ascii_file(path: Path, config: dict[str, Any]) -> Iterator[Spectrum]:
    data = np.loadtxt(
        path,
        delimiter=config.get("delimiter"),
        comments=config.get("comments", "#"),
        skiprows=int(config.get("skiprows", 0)),
    )
    wave = data[:, int(config.get("wavelength_column", 0))]
    flux = data[:, int(config.get("flux_column", 1))]
    raw_uncertainty = None
    if config.get("uncertainty_column") is not None:
        raw_uncertainty = data[:, int(config["uncertainty_column"])]
    yield Spectrum(
        spectrum_id=path.stem,
        wavelength=wave * wavelength_factor(config.get("wavelength_unit", "angstrom")),
        flux=flux * float(config.get("flux_scale", 1.0)),
        uncertainty=_uncertainty(
            raw_uncertainty,
            config.get("uncertainty_kind", "sigma"),
            float(config.get("uncertainty_scale", config.get("flux_scale", 1.0))),
        ),
        redshift=_file_redshift(config),
        metadata={"input_file": str(path)},
    )


def iter_spectrum_files(config: dict[str, Any]) -> Iterator[Spectrum]:
    """Yield one spectrum from each FITS table, FITS image, or ASCII file."""
    file_format = config.get("format", "fits_table")
    readers = {
        "fits_table": _iter_fits_table_file,
        "fits_image": _iter_fits_image_file,
        "ascii": _iter_ascii_file,
    }
    if file_format not in readers:
        raise ConfigError("input.format must be fits_table, fits_image, or ascii")
    reader = readers[file_format]
    for path in _resolve_file_list(config):
        yield from reader(path, config)


def _normalized_miri_segment(value: Any) -> str:
    text = str(value).strip().upper().replace("-", "").replace("_", "")
    for word, letter in (("SHORT", "A"), ("MEDIUM", "B"), ("LONG", "C")):
        text = text.replace(word, letter)
    if len(text) == 2 and text[0] in "1234" and text[1] in "ABC":
        return text
    raise ConfigError("input.segment must be 1A..4C or channel-short/medium/long")


def _select_miri_segment(
    source: dict[str, Any], fit: dict[str, Any]
) -> dict[str, Any]:
    """Select one MIRI band cube that contains the complete fitting window."""
    fits = _fits()
    pattern = source.get("glob")
    if not pattern:
        raise ConfigError("MIRI segmented input requires input.glob or input.path")
    paths = sorted(Path(item) for item in glob.glob(str(pattern)))
    if not paths:
        raise ConfigError(f"MIRI input.glob matched no files: {pattern}")

    lo, hi = map(float, fit["window"])
    if fit.get("frame", "rest") == "rest":
        redshift_factor = 1.0 + float(source["redshift"])
        lo, hi = lo * redshift_factor, hi * redshift_factor
    requested_segment = (
        None if source.get("segment") is None else _normalized_miri_segment(source["segment"])
    )
    candidates: list[tuple[Path, str, float, float]] = []
    coverage_descriptions: list[str] = []
    for path in paths:
        primary = fits.getheader(path, 0)
        science = fits.getheader(path, source.get("flux_hdu", "SCI"))
        channel = str(primary.get("CHANNEL", "")).strip()
        band = str(primary.get("BAND", "")).strip().upper()
        letter = {"SHORT": "A", "MEDIUM": "B", "LONG": "C"}.get(band)
        if channel not in {"1", "2", "3", "4"} or letter is None:
            continue
        segment = f"{channel}{letter}"
        n_pixels = int(science.get("NAXIS3", 0))
        wavelength = linear_wavelength(
            science,
            n_pixels,
            fits_axis=3,
            unit=science.get("CUNIT3", "micron"),
        )
        wave_lo, wave_hi = float(np.min(wavelength)), float(np.max(wavelength))
        coverage_descriptions.append(f"{segment}={wave_lo:.1f}:{wave_hi:.1f} A")
        if requested_segment is not None and segment != requested_segment:
            continue
        if wave_lo <= lo and hi <= wave_hi:
            candidates.append((path, segment, wave_lo, wave_hi))

    if not candidates:
        requested = f" for segment {requested_segment}" if requested_segment else ""
        raise ConfigError(
            f"No MIRI cube{requested} contains observed fitting window "
            f"{lo:.1f}:{hi:.1f} A. Available: " + ", ".join(coverage_descriptions)
        )
    if len(candidates) > 1:
        segments = ", ".join(item[1] for item in candidates)
        raise ConfigError(
            f"Observed fitting window is contained by multiple MIRI segments ({segments}); "
            "set input.segment explicitly"
        )
    selected, segment, wave_lo, wave_hi = candidates[0]
    result = dict(source)
    result["path"] = str(selected)
    result["_segment_selection"] = "explicit" if requested_segment else "wavelength_window"
    result["_segment_coverage_angstrom"] = [wave_lo, wave_hi]
    result["_selected_segment"] = segment
    return result


def iter_spectra(config: dict[str, Any]) -> Iterator[Spectrum]:
    """Dispatch to the configured input adapter."""
    input_config = dict(config["input"])
    input_config["_config_dir"] = config.get("_config_dir", ".")
    if input_config.get("adapter") == "miri" and not input_config.get("path"):
        input_config = _select_miri_segment(input_config, config["fit"])
    kind = input_config["kind"]
    if kind == "cube":
        yield from iter_cube(input_config)
    elif kind == "survey_table":
        yield from iter_survey_table(input_config)
    elif kind == "spectrum_files":
        yield from iter_spectrum_files(input_config)
    else:  # protected by validation
        raise ConfigError(f"Unsupported input kind: {kind}")
