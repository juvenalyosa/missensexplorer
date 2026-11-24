"""MissensExplorer: predict how DNA variants change the encoded protein.

The public API is :func:`predict` (configured with :class:`PredictOptions`) and its shortcut
:func:`predict_text`; see :mod:`missensexplorer.predict`.
"""

# Defined before the import below: api.py and predict.py import it from the package.
__version__ = "2.3.0"

from .predict import PredictOptions, predict, predict_text  # noqa: F401
