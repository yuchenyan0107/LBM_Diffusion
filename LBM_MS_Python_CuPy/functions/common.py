import numpy as np
import matplotlib.pyplot as plt
use_cupy = True

if use_cupy:

    try:
        import cupy as cp
        xp = cp
        print("Using cupy")
    except ImportError:
        cp = None
        xp = np
else:
    cp = None
    xp = np

def to_numpy(arr):
    if cp is not None and hasattr(cp, "asnumpy") and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return np.asarray(arr)

D2Q9_CX = xp.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=xp.float32) # CW: up, ur, r, dr, d.....
D2Q9_CY = xp.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=xp.float32)

w = xp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36], dtype=xp.float32)
OPPOSITE = xp.array([0,5,6,7,8,1,2,3,4], dtype=xp.int32)
theta = 0.5

def safe_divide(numerator, denominator, mask=None):
    """
    Elementwise division that avoids unsupported CuPy `where`/`out` kwargs.
    Returns zero where the mask is False.
    """
    if mask is None:
        mask = denominator != 0
    denom_safe = xp.where(mask, denominator, xp.ones_like(denominator))
    result = numerator / denom_safe
    return xp.where(mask, result, xp.zeros_like(result))