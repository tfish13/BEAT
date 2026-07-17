"""Named instrument adapters layered on BEAT's generic input readers."""

from __future__ import annotations

from typing import Any


class AdapterError(ValueError):
    """Raised when a file does not match its declared instrument adapter."""


SUPPORTED_ADAPTERS = {"miri", "muse", "nirspec"}


NIRSPEC_ETC_LSF = {
    "G235H": {
        "filename": "jwst_nirspec_g235h_disp.fits",
        "sha256": "b699a70a8009f9adc1d1c55516b2731b43c7930de1df10534c64ba56bfb4a413",
        "wavelength_range_micron": [1.5, 3.11001],
    },
    "G395H": {
        "filename": "jwst_nirspec_g395h_disp.fits",
        "sha256": "fc629587007cd8c0e86d803d8f98b5549a2a92993f6f5f73f907db5ceb6b48d8",
        "wavelength_range_micron": [2.5, 6.0],
    },
}


def _nirspec_lsf_metadata(grating: str) -> dict[str, Any]:
    calibration = NIRSPEC_ETC_LSF.get(grating)
    if calibration is None:
        if grating == "PRISM":
            resolving_power = 100.0
        elif grating.endswith("M"):
            resolving_power = 1000.0
        elif grating.endswith("H"):
            resolving_power = 2700.0
        else:
            raise AdapterError(f"Unsupported or missing NIRSpec grating: {grating!r}")
        return {
            "model": "resolving_power",
            "value": resolving_power,
            "source": f"nominal JWST NIRSpec {grating} resolving power",
            "approximation": True,
        }
    return {
        "model": "nirspec_etc_resolving_power",
        "grating": grating,
        "source": "STScI JWST User Documentation ETC dispersion/resolution table",
        "reference_url": (
            "https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/"
            "nirspec-instrumentation/nirspec-dispersers-and-filters"
        ),
        "calibration_file": calibration["filename"],
        "calibration_sha256": calibration["sha256"],
        "wavelength_range_micron": calibration["wavelength_range_micron"],
        "resolution_element_pixels": 2.2,
        "aperture_assumption": "fully illuminated IFU aperture",
        "profile_assumption": "Gaussian with tabulated FWHM=lambda/R",
        "approximation": True,
    }


def apply_adapter_defaults(config: dict[str, Any]) -> None:
    """Apply safe static defaults for a named adapter in-place."""
    source = config.get("input", {})
    adapter = source.get("adapter")
    if adapter is None:
        return
    adapter = str(adapter).lower()
    source["adapter"] = adapter
    if adapter not in SUPPORTED_ADAPTERS:
        raise AdapterError(
            f"Unknown input adapter {adapter!r}; available adapters: "
            + ", ".join(sorted(SUPPORTED_ADAPTERS))
        )

    if adapter == "muse":
        defaults = {
            "kind": "cube",
            "flux_hdu": "DATA",
            "uncertainty_hdu": "STAT",
            "uncertainty_kind": "variance",
            "spectral_axis": 0,
            "wavelength_header_hdu": "DATA",
            "wavelength_unit": "angstrom",
            "flux_scale": 1.0e-20,
            "target_header_hdu": 0,
            "target_header": "OBJECT",
        }
        for key, value in defaults.items():
            source.setdefault(key, value)
        fit = config.setdefault("fit", {})
        fit.setdefault("wavelength_medium", "air")
        # Phase-3 products expose a scalar SPEC_RES value. It is useful as a
        # documented approximation; users can override this with a calibrated
        # wavelength-dependent polynomial or table.
        fit.setdefault("lsf", {"model": "instrument"})
    elif adapter == "nirspec":
        defaults = {
            "kind": "cube",
            "flux_hdu": "SCI",
            "uncertainty_hdu": "ERR",
            "uncertainty_kind": "sigma",
            "spectral_axis": 0,
            "wavelength_header_hdu": "SCI",
            "wavelength_unit": "micron",
            "flux_scale": 1.0,
            "mask_hdu": "DQ",
            # DO_NOT_USE is bit zero. Reading raw storage keeps the large DQ
            # image memory-mapped; its low-order flag bits are unchanged by
            # the FITS unsigned-integer BZERO convention.
            "mask_bits": 1,
            "raw_integer_mask": True,
            "target_header_hdu": 0,
            "target_header": "TARGNAME",
        }
        for key, value in defaults.items():
            source.setdefault(key, value)
        fit = config.setdefault("fit", {})
        fit.setdefault("wavelength_medium", "vacuum")
        fit.setdefault("lsf", {"model": "instrument"})
    elif adapter == "miri":
        defaults = {
            "kind": "cube",
            "flux_hdu": "SCI",
            "uncertainty_hdu": "ERR",
            "uncertainty_kind": "sigma",
            "spectral_axis": 0,
            "wavelength_header_hdu": "SCI",
            "wavelength_unit": "micron",
            "flux_scale": 1.0,
            "mask_hdu": "DQ",
            "mask_bits": 1,
            "raw_integer_mask": True,
            "target_header_hdu": 0,
            "target_header": "TARGNAME",
        }
        for key, value in defaults.items():
            source.setdefault(key, value)
        fit = config.setdefault("fit", {})
        fit.setdefault("wavelength_medium", "vacuum")
        fit.setdefault("lsf", {"model": "instrument"})


def _find_lsf_calibration(header: Any) -> dict[str, Any]:
    for key, value in header.items():
        if str(value).strip().upper() != "LSF_PROFILE" or not key.endswith(" CATG"):
            continue
        prefix = key[: -len(" CATG")]
        return {
            "lsf_calibration_name": header.get(prefix + " NAME"),
            "lsf_calibration_md5": header.get(prefix + " DATAMD5"),
        }
    return {}


def inspect_muse_hdul(hdul: Any, source: dict[str, Any]) -> dict[str, Any]:
    """Validate an ESO Phase-3 MUSE cube and return provenance metadata."""
    primary = hdul[0].header
    instrument = str(primary.get("INSTRUME", "")).strip().upper()
    if instrument != "MUSE":
        raise AdapterError(
            f"MUSE adapter expected primary INSTRUME='MUSE', found {instrument!r}"
        )
    product = str(primary.get("PRODCATG", "")).strip().upper()
    if product and product != "SCIENCE.CUBE.IFS":
        raise AdapterError(
            f"MUSE adapter expected PRODCATG='SCIENCE.CUBE.IFS', found {product!r}"
        )

    try:
        data_hdu = hdul[source["flux_hdu"]]
        stat_hdu = hdul[source["uncertainty_hdu"]]
    except (KeyError, IndexError) as exc:
        raise AdapterError("MUSE cube must contain DATA and STAT extensions") from exc
    if data_hdu.data is None or data_hdu.data.ndim != 3:
        raise AdapterError("MUSE DATA extension must contain a three-dimensional cube")
    if stat_hdu.data is None or stat_hdu.data.shape != data_hdu.data.shape:
        raise AdapterError("MUSE STAT extension must match the DATA cube shape")
    ctype = str(data_hdu.header.get("CTYPE3", "")).upper()
    if ctype != "AWAV":
        raise AdapterError(
            f"MUSE adapter requires an air-wavelength AWAV axis, found CTYPE3={ctype!r}"
        )

    pipeline_versions = []
    for key, value in primary.items():
        if key.startswith("ESO PRO REC") and key.endswith(" PIPE ID"):
            pipeline_versions.append(str(value))
    resolving_power = primary.get("SPEC_RES")
    if resolving_power is None or float(resolving_power) <= 0:
        raise AdapterError("MUSE Phase-3 primary header is missing positive SPEC_RES")

    metadata: dict[str, Any] = {
        "adapter": "muse",
        "instrument": "MUSE",
        "product_category": product or None,
        "pipeline_versions": sorted(set(pipeline_versions)),
        "wavelength_medium": "air",
        "flux_unit_raw": data_hdu.header.get("BUNIT", primary.get("BUNIT")),
        "uncertainty_unit_raw": stat_hdu.header.get("BUNIT"),
        "spectral_sampling_angstrom": float(
            data_hdu.header.get("CDELT3", data_hdu.header.get("CD3_3"))
        ),
        "instrument_lsf": {
            "model": "resolving_power",
            "value": float(resolving_power),
            "source": "MUSE Phase-3 primary header SPEC_RES",
            "approximation": True,
        },
    }
    metadata.update(_find_lsf_calibration(primary))
    return metadata


def inspect_nirspec_hdul(hdul: Any, source: dict[str, Any]) -> dict[str, Any]:
    """Validate a JWST NIRSpec Stage-3 IFU cube and return provenance."""
    primary = hdul[0].header
    if str(primary.get("TELESCOP", "")).strip().upper() != "JWST":
        raise AdapterError("NIRSpec adapter requires primary TELESCOP='JWST'")
    if str(primary.get("INSTRUME", "")).strip().upper() != "NIRSPEC":
        raise AdapterError("NIRSpec adapter requires primary INSTRUME='NIRSPEC'")
    if str(primary.get("EXP_TYPE", "")).strip().upper() != "NRS_IFU":
        raise AdapterError("NIRSpec cube must have EXP_TYPE='NRS_IFU'")
    data_model = str(primary.get("DATAMODL", "")).strip()
    if data_model and data_model != "IFUCubeModel":
        raise AdapterError(
            f"NIRSpec adapter expected DATAMODL='IFUCubeModel', found {data_model!r}"
        )

    try:
        science = hdul[source["flux_hdu"]]
        error = hdul[source["uncertainty_hdu"]]
        dq = hdul[source["mask_hdu"]]
    except (KeyError, IndexError) as exc:
        raise AdapterError("NIRSpec s3d cube must contain SCI, ERR, and DQ extensions") from exc
    expected_shape = (
        science.header.get("NAXIS3"),
        science.header.get("NAXIS2"),
        science.header.get("NAXIS1"),
    )
    if science.header.get("NAXIS") != 3 or any(value is None for value in expected_shape):
        raise AdapterError("NIRSpec SCI extension must be a three-dimensional cube")
    for extension in (error, dq):
        shape = (
            extension.header.get("NAXIS3"),
            extension.header.get("NAXIS2"),
            extension.header.get("NAXIS1"),
        )
        if shape != expected_shape:
            raise AdapterError("NIRSpec ERR and DQ extensions must match SCI")
    if str(science.header.get("CTYPE3", "")).strip().upper() != "WAVE":
        raise AdapterError("NIRSpec SCI extension requires CTYPE3='WAVE'")
    if str(science.header.get("CUNIT3", "")).strip().lower() not in {"um", "micron"}:
        raise AdapterError("NIRSpec SCI wavelength axis must be calibrated in microns")
    if str(science.header.get("BUNIT", "")).strip().lower() != "mjy/sr":
        raise AdapterError("NIRSpec SCI extension must use BUNIT='MJy/sr'")

    grating = str(primary.get("GRATING", "")).strip().upper()
    instrument_lsf = _nirspec_lsf_metadata(grating)
    sampling = science.header.get("CDELT3")
    metadata: dict[str, Any] = {
        "adapter": "nirspec",
        "telescope": "JWST",
        "instrument": "NIRSpec",
        "exposure_type": primary.get("EXP_TYPE"),
        "data_model": data_model or None,
        "pipeline_version": primary.get("CAL_VER"),
        "crds_context": primary.get("CRDS_CTX"),
        "grating": grating,
        "filter": primary.get("FILTER"),
        "source_type": science.header.get("SRCTYPE"),
        "wavelength_medium": "vacuum",
        "flux_unit_raw": science.header.get("BUNIT"),
        "uncertainty_unit_raw": error.header.get("BUNIT"),
        "pixel_area_steradian": science.header.get("PIXAR_SR"),
        "pixel_area_arcsec2": science.header.get("PIXAR_A2"),
        "spectral_sampling_angstrom": (
            None if sampling is None else float(sampling) * 1.0e4
        ),
        "wavelength_table_used": "WCS-TABLE" in hdul,
        "instrument_lsf": instrument_lsf,
        "disperser_reference": primary.get("R_DISPER"),
        "cube_parameters_reference": primary.get("R_CUBPAR"),
        "photom_reference": primary.get("R_PHOTOM"),
    }
    return metadata


def inspect_miri_hdul(hdul: Any, source: dict[str, Any]) -> dict[str, Any]:
    """Validate a JWST MIRI MRS band cube and return segment provenance."""
    primary = hdul[0].header
    if str(primary.get("TELESCOP", "")).strip().upper() != "JWST":
        raise AdapterError("MIRI adapter requires primary TELESCOP='JWST'")
    if str(primary.get("INSTRUME", "")).strip().upper() != "MIRI":
        raise AdapterError("MIRI adapter requires primary INSTRUME='MIRI'")
    if str(primary.get("EXP_TYPE", "")).strip().upper() != "MIR_MRS":
        raise AdapterError("MIRI adapter requires EXP_TYPE='MIR_MRS'")
    data_model = str(primary.get("DATAMODL", "")).strip()
    if data_model and data_model != "IFUCubeModel":
        raise AdapterError(
            f"MIRI adapter expected DATAMODL='IFUCubeModel', found {data_model!r}"
        )
    try:
        science = hdul[source["flux_hdu"]]
        error = hdul[source["uncertainty_hdu"]]
        dq = hdul[source["mask_hdu"]]
    except (KeyError, IndexError) as exc:
        raise AdapterError("MIRI s3d cube must contain SCI, ERR, and DQ extensions") from exc
    expected_shape = tuple(
        science.header.get(f"NAXIS{axis}") for axis in (3, 2, 1)
    )
    if science.header.get("NAXIS") != 3 or any(value is None for value in expected_shape):
        raise AdapterError("MIRI SCI extension must be a three-dimensional cube")
    for extension in (error, dq):
        shape = tuple(extension.header.get(f"NAXIS{axis}") for axis in (3, 2, 1))
        if shape != expected_shape:
            raise AdapterError("MIRI ERR and DQ extensions must match SCI")
    if str(science.header.get("CTYPE3", "")).strip().upper() != "WAVE":
        raise AdapterError("MIRI SCI extension requires CTYPE3='WAVE'")
    if str(science.header.get("CUNIT3", "")).strip().lower() not in {"um", "micron"}:
        raise AdapterError("MIRI SCI wavelength axis must be calibrated in microns")
    if str(science.header.get("BUNIT", "")).strip().lower() != "mjy/sr":
        raise AdapterError("MIRI SCI extension must use BUNIT='MJy/sr'")
    channel = str(primary.get("CHANNEL", "")).strip()
    band = str(primary.get("BAND", "")).strip().upper()
    if channel not in {"1", "2", "3", "4"} or band not in {
        "SHORT", "MEDIUM", "LONG"
    }:
        raise AdapterError(f"Invalid MIRI MRS channel/band: {channel!r}/{band!r}")
    band_letter = {"SHORT": "A", "MEDIUM": "B", "LONG": "C"}[band]
    sampling = science.header.get("CDELT3")
    return {
        "adapter": "miri",
        "telescope": "JWST",
        "instrument": "MIRI",
        "exposure_type": primary.get("EXP_TYPE"),
        "data_model": data_model or None,
        "pipeline_version": primary.get("CAL_VER"),
        "crds_context": primary.get("CRDS_CTX"),
        "channel": int(channel),
        "band": band,
        "segment": f"{channel}{band_letter}",
        "segment_selection": source.get("_segment_selection", "path"),
        "segment_coverage_angstrom": source.get("_segment_coverage_angstrom"),
        "source_type": science.header.get("SRCTYPE"),
        "wavelength_medium": "vacuum",
        "flux_unit_raw": science.header.get("BUNIT"),
        "uncertainty_unit_raw": error.header.get("BUNIT"),
        "pixel_area_steradian": science.header.get("PIXAR_SR"),
        "pixel_area_arcsec2": science.header.get("PIXAR_A2"),
        "spectral_sampling_angstrom": (
            None if sampling is None else float(sampling) * 1.0e4
        ),
        "wavelength_table_used": "WCS-TABLE" in hdul,
        "instrument_lsf": {
            "model": "polynomial_resolving_power",
            # STScI's approximate commissioning fit, with wavelength in micron.
            "coefficients": [4603.0, -128.0],
            "reference_angstrom": 0.0,
            "scale_angstrom": 1.0e4,
            "source": "STScI approximate MIRI MRS in-flight R(lambda)",
            "approximation": True,
        },
        "cube_parameters_reference": primary.get("R_CUBPAR"),
        "distortion_reference": primary.get("R_DISTOR"),
        "fringe_reference": primary.get("R_FRINGE"),
        "photom_reference": primary.get("R_PHOTOM"),
    }


def inspect_adapter_hdul(hdul: Any, source: dict[str, Any]) -> dict[str, Any]:
    adapter = source.get("adapter")
    if adapter is None:
        return {}
    if adapter == "muse":
        return inspect_muse_hdul(hdul, source)
    if adapter == "nirspec":
        return inspect_nirspec_hdul(hdul, source)
    if adapter == "miri":
        return inspect_miri_hdul(hdul, source)
    raise AdapterError(f"Unsupported adapter: {adapter}")
