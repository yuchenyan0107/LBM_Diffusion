import numpy as np
import matplotlib.pyplot as plt
use_GPU = True
use_GPU = False

cp = None
mx = None
xp = np

if use_GPU:
    try:
        import cupy as cp
        xp = cp
        print("Using cupy")
    except ImportError:
        cp = None
        try:
            import mlx.core as mx
            xp = mx
            print("Using mlx")
        except ImportError:
            mx = None
            xp = np
            print("using numpy")
else:
    print("Using numpy")
DTYPE = xp.float32


def to_numpy(arr, *, copy=True):
    """
    Convert CuPy / MLX / NumPy-ish arrays to a NumPy ndarray.

    Parameters
    ----------
    copy : bool
        - For MLX: if False, tries to create a NumPy *view* (zero-copy).
        - For CuPy: conversion always copies from device to host.
    """
    # --- CuPy -> NumPy ---
    if cp is not None and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)

    # --- MLX -> NumPy ---
    if mx is not None:
        mlx_array_type = getattr(mx, "array", None)  # mlx.core.array type (callable/class)
        if mlx_array_type is not None and isinstance(arr, mlx_array_type):
            # MLX docs: convert via np.array(...); copy=False yields a view when possible
            return np.array(arr, copy=copy)

    # --- Everything else ---
    return np.asarray(arr)


D2Q9_CX = xp.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=xp.float32) # CW: up, ur, r, dr, d.....
D2Q9_CY = xp.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=xp.float32)

w = xp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36], dtype=xp.float32)
OPPOSITE = xp.array([0,5,6,7,8,1,2,3,4], dtype=xp.float32)
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