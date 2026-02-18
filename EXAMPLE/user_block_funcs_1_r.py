import numpy as np

####################################################################################################
####################################################################################################
####################### define the blocked entry-evaluation function #######################

def compute_block(rows, cols, meta):
    xyz = meta["coordinates"]

    R = xyz[rows]          # shape (nr, D)
    C = xyz[cols]          # shape (nc, D)

    # Pairwise differences: (nr, nc, D)
    diff = R[:, None, :] - C[None, :, :]

    # Squared distances: (nr, nc)
    dist2 = np.einsum("ijk,ijk->ij", diff, diff)

    # Avoid division by zero on diagonal
    out = np.empty(dist2.shape, dtype=np.float64, order="F")
    mask = dist2 > 0
    out[mask] = 1.0 / np.sqrt(dist2[mask])
    out[~mask] = 1.0

    return out

