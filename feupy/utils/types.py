# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for types validation."""

from typing import Annotated, List, Union
from pydantic.functional_validators import BeforeValidator
from feupy.analysis.irfs import Irfs


__all__ = [
    "IrfType",
]


def validate_irf(v):
    if v not in Irfs.IRFS_OPTIONS:
        raise ValueError(f"Invalid IRF option: {v!r}. Choose from: {Irfs.IRFS_OPTIONS!r}")
    return v
    
IrfType = Annotated[
    Union[str, List[str]],
    BeforeValidator(validate_irf),
]
