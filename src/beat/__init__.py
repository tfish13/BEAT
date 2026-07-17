"""BEAT: Bayesian model selection for astronomical emission lines."""

from .config import ConfigError, load_config, validate_config
from .spectrum import Spectrum

__all__ = ["ConfigError", "Spectrum", "load_config", "validate_config"]
__version__ = "2.0.0a1"
