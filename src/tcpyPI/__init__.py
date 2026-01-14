"""
tcpyPI: Tropical Cyclone Potential Intensity Calculations in Python
----------------------------
tcpyPI, 'pyPI' for short, is a set of functions, scripts, and notebooks that compute and validate tropical cyclone potential intensity (PI) calculations in Python. It is a port of the Bister and Emanuel (2002) algorithm (Bister and Emanuel 2002) which was originally written in FORTRAN---and then converted to MATLAB---by Prof. Kerry Emanuel (MIT).
"""

__author__ = "Daniel M. Gilford <dgilford@climatecentral.org>"

# TODO: Re-factor module-wide implicit imports
from .pi import *
from .constants import *
from .utilities import *

# User-facing analyses API
from .analyses import log_decompose_pi, pi_log_decomposition

# Ventilation Index (TE12) diagnostics
from .vi import ventilation_index, vi
from .vi import entropy_deficit_te12, entropy_deficit_te12_from_profile
