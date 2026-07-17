"""Dependency-backed smoke tests for FITS adapters and the real sampler."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

try:
    from astropy.io import fits
    import ultranest  # noqa: F401 - confirms the real sampling dependency exists

    HAS_FITS_STACK = True
except ImportError:
    HAS_FITS_STACK = False

from beat.data import iter_cube, iter_spectra, iter_survey_table
from beat.fitting import fit_spectrum
from beat.model import gaussian_integrated


@unittest.skipUnless(HAS_FITS_STACK, "Astropy and UltraNest are required")
class FitsAndUltraNestTests(unittest.TestCase):
    def test_miri_glob_selects_segment_containing_fit_window(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)

            def write_segment(channel, band, start_micron, stop_micron):
                wave = np.linspace(start_micron, stop_micron, 101)
                cube = np.ones((wave.size, 3, 4), dtype=np.float32)
                primary = fits.PrimaryHDU()
                primary.header["TELESCOP"] = "JWST"
                primary.header["INSTRUME"] = "MIRI"
                primary.header["EXP_TYPE"] = "MIR_MRS"
                primary.header["DATAMODL"] = "IFUCubeModel"
                primary.header["TARGNAME"] = "TEST MIRI"
                primary.header["CHANNEL"] = str(channel)
                primary.header["BAND"] = band
                primary.header["CAL_VER"] = "2.0.1"
                science = fits.ImageHDU(cube, name="SCI")
                science.header["CTYPE3"] = "WAVE"
                science.header["CUNIT3"] = "um"
                science.header["CRVAL3"] = wave[0]
                science.header["CRPIX3"] = 1.0
                science.header["CDELT3"] = wave[1] - wave[0]
                science.header["BUNIT"] = "MJy/sr"
                science.header["PIXAR_SR"] = 9.4e-13
                error = fits.ImageHDU(np.full_like(cube, 0.5), name="ERR")
                error.header["BUNIT"] = "MJy/sr"
                dq = fits.ImageHDU(np.zeros(cube.shape, dtype=np.uint32), name="DQ")
                path = root / f"target_ch{channel}-{band.lower()}_s3d.fits"
                fits.HDUList([primary, science, error, dq]).writeto(path)

            write_segment(3, "MEDIUM", 13.3, 15.6)
            write_segment(3, "LONG", 15.4, 18.0)
            config = {
                "input": {
                    "adapter": "miri",
                    "kind": "cube",
                    "glob": str(root / "target_ch*.fits"),
                    "flux_hdu": "SCI",
                    "uncertainty_hdu": "ERR",
                    "uncertainty_kind": "sigma",
                    "mask_hdu": "DQ",
                    "mask_bits": 1,
                    "raw_integer_mask": True,
                    "wavelength_header_hdu": "SCI",
                    "wavelength_unit": "micron",
                    "redshift": 0.0,
                    "x_range": [1, 2],
                    "y_range": [1, 2],
                    "target_header": "TARGNAME",
                },
                "fit": {"frame": "rest", "window": [140000.0, 145000.0]},
            }
            spectra = list(iter_spectra(config))
            self.assertEqual(len(spectra), 1)
            self.assertEqual(spectra[0].metadata["segment"], "3B")
            self.assertEqual(
                spectra[0].metadata["segment_selection"], "wavelength_window"
            )
            self.assertIn("medium", spectra[0].metadata["input_file"])

    def test_nirspec_adapter_reads_wavelength_table_and_dq_lazily(self) -> None:
        wave_micron = np.linspace(1.70, 1.80, 31)
        cube = np.ones((wave_micron.size, 3, 4), dtype=np.float32)
        dq_cube = np.zeros(cube.shape, dtype=np.uint32)
        dq_cube[5, 1, 2] = 1
        with TemporaryDirectory() as directory:
            path = Path(directory) / "nirspec_s3d.fits"
            primary = fits.PrimaryHDU()
            primary.header["TELESCOP"] = "JWST"
            primary.header["INSTRUME"] = "NIRSPEC"
            primary.header["EXP_TYPE"] = "NRS_IFU"
            primary.header["DATAMODL"] = "IFUCubeModel"
            primary.header["TARGNAME"] = "TEST NIRSPEC"
            primary.header["GRATING"] = "G235H"
            primary.header["FILTER"] = "F170LP"
            primary.header["CAL_VER"] = "2.0.1"
            primary.header["R_DISPER"] = "crds://test_disperser.asdf"
            science = fits.ImageHDU(cube, name="SCI")
            science.header["CTYPE3"] = "WAVE"
            science.header["CUNIT3"] = "um"
            science.header["CRVAL3"] = wave_micron[0]
            science.header["CRPIX3"] = 1.0
            science.header["CDELT3"] = wave_micron[1] - wave_micron[0]
            science.header["BUNIT"] = "MJy/sr"
            science.header["PIXAR_SR"] = 2.35e-13
            error = fits.ImageHDU(np.full_like(cube, 0.25), name="ERR")
            error.header["BUNIT"] = "MJy/sr"
            dq = fits.ImageHDU(dq_cube, name="DQ")
            wave_column = fits.Column(
                name="wavelength",
                format=f"{wave_micron.size}D",
                array=[wave_micron],
            )
            wave_table = fits.BinTableHDU.from_columns(
                [wave_column], name="WCS-TABLE"
            )
            fits.HDUList([primary, science, error, dq, wave_table]).writeto(path)

            spectra = list(
                iter_cube(
                    {
                        "adapter": "nirspec",
                        "path": str(path),
                        "flux_hdu": "SCI",
                        "uncertainty_hdu": "ERR",
                        "uncertainty_kind": "sigma",
                        "mask_hdu": "DQ",
                        "mask_bits": 1,
                        "raw_integer_mask": True,
                        "wavelength_header_hdu": "SCI",
                        "wavelength_unit": "micron",
                        "redshift": 0.01,
                        "x_range": [2, 3],
                        "y_range": [1, 2],
                        "target_header": "TARGNAME",
                        "target_header_hdu": 0,
                    }
                )
            )
            self.assertEqual(len(spectra), 1)
            spectrum = spectra[0]
            np.testing.assert_allclose(spectrum.wavelength, wave_micron * 1.0e4)
            self.assertTrue(spectrum.mask[5])
            self.assertEqual(spectrum.metadata["grating"], "G235H")
            self.assertEqual(
                spectrum.metadata["instrument_lsf"]["model"],
                "nirspec_etc_resolving_power",
            )
            self.assertEqual(
                spectrum.metadata["instrument_lsf"]["calibration_file"],
                "jwst_nirspec_g235h_disp.fits",
            )
            self.assertTrue(spectrum.metadata["wavelength_table_used"])

    def test_muse_adapter_validates_and_preserves_provenance(self) -> None:
        wave = np.linspace(4800.0, 5100.0, 61)
        cube = np.ones((wave.size, 4, 5), dtype=np.float32)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "muse.fits"
            primary = fits.PrimaryHDU()
            primary.header["INSTRUME"] = "MUSE"
            primary.header["PRODCATG"] = "SCIENCE.CUBE.IFS"
            primary.header["OBJECT"] = "TEST GALAXY"
            primary.header["SPEC_RES"] = 3000.0
            primary.header["HIERARCH ESO PRO REC1 PIPE ID"] = "muse/2.8"
            primary.header["HIERARCH ESO PRO REC1 CAL1 CATG"] = "LSF_PROFILE"
            primary.header["HIERARCH ESO PRO REC1 CAL1 NAME"] = "test_lsf.fits"
            primary.header["HIERARCH ESO PRO REC1 CAL1 DATAMD5"] = "abc123"
            data = fits.ImageHDU(cube, name="DATA")
            data.header["CTYPE3"] = "AWAV"
            data.header["CRVAL3"] = wave[0]
            data.header["CRPIX3"] = 1.0
            data.header["CDELT3"] = wave[1] - wave[0]
            data.header["BUNIT"] = "1e-20 erg/s/cm2/Angstrom"
            stat = fits.ImageHDU(np.full_like(cube, 4.0), name="STAT")
            fits.HDUList([primary, data, stat]).writeto(path)

            spectra = list(
                iter_cube(
                    {
                        "adapter": "muse",
                        "path": str(path),
                        "flux_hdu": "DATA",
                        "uncertainty_hdu": "STAT",
                        "uncertainty_kind": "variance",
                        "wavelength_header_hdu": "DATA",
                        "flux_scale": 1.0e-20,
                        "redshift": 0.01,
                        "x_range": [1, 3],
                        "y_range": [2, 4],
                        "target_header": "OBJECT",
                        "target_header_hdu": 0,
                    }
                )
            )
            self.assertEqual(len(spectra), 4)
            self.assertTrue(spectra[0].spectrum_id.startswith("TEST GALAXY_"))
            self.assertEqual(spectra[0].metadata["instrument"], "MUSE")
            self.assertTrue(spectra[0].metadata["instrument_lsf"]["approximation"])
            self.assertEqual(spectra[0].metadata["lsf_calibration_name"], "test_lsf.fits")
            self.assertAlmostEqual(spectra[0].uncertainty[0], 2.0e-20)

    def test_cube_and_survey_table_with_real_nested_sampling(self) -> None:
        np.random.seed(42)
        rng = np.random.default_rng(42)
        wave = np.linspace(4980.0, 5030.0, 101)
        noise = 0.15
        injected_flux = 20.0
        flux = gaussian_integrated(wave, 5006.84, 1.2, injected_flux)
        flux += rng.normal(0.0, noise, wave.size)

        with TemporaryDirectory() as directory:
            root = Path(directory)
            cube_path = root / "cube.fits"
            cube_flux = np.broadcast_to(flux[:, None, None], (wave.size, 2, 2)).copy()
            cube_hdu = fits.ImageHDU(cube_flux)
            cube_hdu.header["CRVAL3"] = wave[0]
            cube_hdu.header["CRPIX3"] = 1.0
            cube_hdu.header["CDELT3"] = wave[1] - wave[0]
            fits.HDUList(
                [fits.PrimaryHDU(), cube_hdu, fits.ImageHDU(np.full_like(cube_flux, noise))]
            ).writeto(cube_path)
            cube_spectra = list(
                iter_cube(
                    {
                        "path": str(cube_path),
                        "flux_hdu": 1,
                        "uncertainty_hdu": 2,
                        "uncertainty_kind": "sigma",
                        "redshift": 0.0,
                    }
                )
            )
            self.assertEqual(len(cube_spectra), 4)
            np.testing.assert_allclose(cube_spectra[0].wavelength, wave)

            survey_path = root / "survey.fits"
            columns = [
                fits.Column(name="TARGETID", format="16A", array=["target_a", "target_b"]),
                fits.Column(name="Z", format="D", array=[0.0, 0.0]),
                fits.Column(name="WAVE", format=f"{wave.size}D", array=[wave, wave]),
                fits.Column(name="FLUX", format=f"{wave.size}D", array=[flux, flux]),
                fits.Column(
                    name="IVAR",
                    format=f"{wave.size}D",
                    array=[np.full(wave.size, noise**-2)] * 2,
                ),
            ]
            fits.HDUList(
                [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="SPECTRA")]
            ).writeto(survey_path)
            survey_spectra = list(
                iter_survey_table(
                    {
                        "path": str(survey_path),
                        "hdu": "SPECTRA",
                        "id_column": "TARGETID",
                        "redshift_column": "Z",
                        "wavelength_column": "WAVE",
                        "flux_column": "FLUX",
                        "uncertainty_column": "IVAR",
                        "uncertainty_kind": "inverse_variance",
                    }
                )
            )
            self.assertEqual([item.spectrum_id for item in survey_spectra], ["target_a", "target_b"])
            np.testing.assert_allclose(survey_spectra[0].uncertainty, noise)

            fit = {
                "frame": "rest",
                "window": [4980.0, 5030.0],
                "continuum": {"degree": 0, "windows": [[4980.0, 4990.0], [5020.0, 5030.0]]},
                "kinematics": {
                    "max_components": 1,
                    "velocity_kms": [-250.0, 250.0],
                    "sigma_kms": [25.0, 200.0],
                },
                "lines": [{"name": "oiii5007", "wavelength": 5006.84}],
                "selection": {"delta_logz": 2.0},
                "sampling": {
                    "min_num_live_points": 40,
                    "min_ess": 40,
                    "dlogz": 5.0,
                    "show_status": False,
                    "stepsampler": "none",
                },
                "minimum_valid_pixels": 20,
            }
            # UltraNest writes progress messages even with show_status=False.
            # Keep the regression suite readable while still using the real sampler.
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                result = fit_spectrum(survey_spectra[0], fit)
            self.assertEqual(result["selected_components"], 1)
            self.assertIn(result["selection_status"], {"accepted", "ambiguous"})
            self.assertIn("selection_diagnostics", result)
            self.assertFalse(result["selection_audit"]["performed"])
            recovered = result["components"][0]["lines"]["oiii5007"]["flux"]
            self.assertAlmostEqual(recovered, injected_flux, delta=3.0)


if __name__ == "__main__":
    unittest.main()
