import numpy as np


def naninterp(field):

    coefficients = [np.finfo(float).eps, np.finfo(float).eps, 1/2, 4/3, 1/2]

    n_miss = len(np.where(np.isnan(field))[0])
    szM = field.shape
    szMf2 = szM[1]+4
    isvec = np.array([1 if x == 0 else 0 for x in szM])
    v4 = [x for x in range(-1, 3)]
    o4 = np.ones((1, 4))
    o5 = np.ones((5, 1))
    o2 = np.ones((2, 1))
    omiss = np.ones((n_miss, 1))
    wsum = np.zeros((n_miss, 2))
    a = wsum.copy()

    for jj in np.where(isvec == 0)[0]:
        bM = np.zeros((2, szM[2-jj-1]))
        Mf = field
        if jj == 1:
            Mf = Mf.T
        Mf = np.vstack((bM, np.isnan(Mf), bM))
        Mf = Mf.flatten('F')
        miss = np.where(Mf == 1)[0]
        exis = np.where(Mf == 0)
        Mf = np.cumsum(np.logical_not(Mf))
        i_m = Mf[miss]

        Mf = field
        if jj == 1:
            Mf = Mf.T
        bMnan = bM.copy()
        bMnan[bM == 0] = np.nan
        Mf = np.vstack((bMnan, Mf, bMnan))
        I = np.tile(i_m[:, np.newaxis], (1, 4))+np.tile(v4, (n_miss, 1))
        I = exis[0][np.unravel_index(I - 1, exis[0].shape)]

        W = np.tile(miss[:, np.newaxis], (1, 4)) - I
        A = np.zeros((n_miss, 5))
        # Mf[np.unravel_index(I[:,1:3],Mf.shape,'F')]
        col0 = np.unravel_index(I[:, 0], Mf.shape, 'F')
        col0 = Mf[col0[0], col0[1]]
        col1 = np.unravel_index(I[:, 1], Mf.shape, 'F')
        col1 = Mf[col1[0], col1[1]]
        col2 = np.unravel_index(I[:, 2], Mf.shape, 'F')
        col2 = Mf[col2[0], col2[1]]
        col3 = np.unravel_index(I[:, 3], Mf.shape, 'F')
        col3 = Mf[col3[0], col3[1]]
        A[:, 0:2] = np.reshape(np.hstack((col1, col2)), (n_miss, 2), 'F')
        A[:, 2:5] = W[:, 0:3]-W[:, 1:4]
        A[:, 2:5] = np.multiply(np.reshape(np.hstack((col1, col2, col3)), (n_miss, 3), 'F'),W[:, 0:3])
        A[:, 2:5] = A[:, 2:5] - np.multiply(np.vstack((col0, col1, col2)).T, W[:, 1:4])
        A[:, 2:5] = np.divide(A[:, 2:5], (W[:, 0:3]-W[:, 1:4]))

        W = np.hstack((np.abs(W[:, 1:3]), np.abs(W[:, 0:3])+np.abs(W[:, 1:4])))
        W = np.divide(np.logical_not(np.isnan(A)), W)
        W = np.multiply(W, np.tile(coefficients, (n_miss, 1)))
        wsum[:, jj] = np.sum(W, 1)
        wsum[:, jj] = wsum[:, jj] - (wsum[:, jj] == 0)
        W = np.divide(W, np.tile(wsum[:, jj], (5, 1)).T)
        A[np.isnan(A)] = 0
        A = np.multiply(A, W)
        a[:, jj] = (A@o5).T

    exis = np.ceil((miss+1)/szMf2)
    exis = exis + (miss+1 - (exis - 1) * szMf2) * szMf2
    i_m = np.argsort(exis)
    exis = np.sort(exis)
    wsum[:, 0] = wsum[i_m, 0]
    asort = a.copy()
    asort[i_m, 0] = a[:, 0]

    wsum = wsum + (wsum == -1)
    exis = wsum@o2
    i_m = exis == 0
    exis[i_m] = exis[i_m] + np.nan
    wsum = np.divide(wsum, np.tile(exis, (1, 2)))
    exis = np.multiply(asort, wsum)@o2

    res = np.unravel_index(miss, Mf.shape, 'F')
    Mf[res[0], res[1]] = exis.T
    Mf = Mf[2:-2, :]
    if jj == 1:
        Mf = Mf.T

    return Mf
