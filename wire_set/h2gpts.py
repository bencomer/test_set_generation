import numpy as np

def h2gpts(h, cell_cv, idiv=4):
    """Convert grid spacing to number of grid points divisible by idiv.
    Taken from GPAW:
        https://gitlab.com/gpaw/gpaw/blob/master/gpaw/utilities/__init__.py

    Note that units of h and cell_cv must match!

    h: float
        Desired grid spacing in.
    cell_cv: 3x3 ndarray
        Unit cell.
    """

    L_c = (np.linalg.inv(cell_cv)**2).sum(0)**-0.5
    print(L_c)
    return np.maximum(idiv, (L_c / h / idiv + 0.5).astype(int) * idiv)


cell = np.array([[1.5119999999999998, 5.642860821044143, 0.0], [0.0, 40.0, 0.0], [0.0, 0.0, 15.0]])


grid = h2gpts(0.12, cell, idiv=1)
print(grid)
