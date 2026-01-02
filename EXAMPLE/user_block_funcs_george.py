import numpy as np

def compute_block(rows, cols, meta):
    xyz = meta["coordinates"]
    K = meta["kernel"]
    yerr = meta["yerr"]

    R = xyz[rows]          # shape (nr, D)
    C = xyz[cols]          # shape (nc, D)

    out = K.kernel.value_general(R,C).astype(np.float64,order="F")

    col_pos = {c: j for j, c in enumerate(cols)}
    for i, ri in enumerate(rows):
        j = col_pos.get(ri)
        if j is not None:
            out[i, j] += yerr[ri] ** 2

    return out





