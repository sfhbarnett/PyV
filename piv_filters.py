import numpy as np


def localfilt(x, y, u, v, threshold):
    m = 3
    nu = np.empty((int(u.shape[0]+2*np.floor(m/2)), int(u.shape[1]+2*np.floor(m/2))))
    nv = np.empty((int(u.shape[0]+2*np.floor(m/2)), int(u.shape[1]+2*np.floor(m/2))))
    nu.fill(np.nan)
    nv.fill(np.nan)
    nu[int(np.floor(m/2)):-int(np.floor(m/2)), int(np.floor(m/2)):-int(np.floor(m/2))] = u
    nv[int(np.floor(m/2)):-int(np.floor(m/2)), int(np.floor(m/2)):-int(np.floor(m/2))] = v

    U2 = nu + nv * 1j
    histostd = np.zeros(nu.shape, dtype='complex')
    histo = np.zeros(nu.shape, dtype='complex')
    ma, na = U2.shape
    INx = np.zeros(nu.shape)
    for ii in range(m-2, na-m+2, 1):
        for jj in range(m-2, ma-m+2, 1):
            if INx[jj, ii] != 1:
                tmp = U2[jj - int(np.floor(m / 2)):jj + int(np.floor(m / 2)) + 1,
                      ii - int(np.floor(m / 2)):ii + int(np.floor(m / 2)) + 1].copy()
                tmp[int(np.floor(m/2)), int(np.floor(m/2))] = np.nan+np.nan*1j
                notnans = np.invert(np.isnan(tmp))
                if np.any(notnans):
                    t = tmp.flatten()
                    tabs = np.abs(t)
                    if np.sum(notnans) % 2 == 1:
                        histo[jj, ii] = t[np.where(tabs == np.nanmedian(tabs))][0]
                        histostd[jj, ii] = np.std(np.real(tmp[notnans]))+np.std(np.imag(tmp[notnans]))*1j
                    else:
                        notnans = np.invert(np.isnan(t))
                        complex = t[notnans]
                        idx = np.argsort(tabs[notnans])
                        complex = complex[idx]
                        middle = int(len(complex)/2)
                        left = complex[middle-1]
                        right = complex[middle]
                        numberreal = np.mean([np.real(left), np.real(right)])
                        numberimag = np.mean([np.imag(left), np.imag(right)])
                        histo[jj, ii] = numberreal+numberimag*1j
                        histostd[jj, ii] = np.std(np.real(complex))+np.std(np.imag(complex))*1j
                else:
                    histo[jj, ii] = np.nan+np.nan*1j
                    histostd[jj, ii] = np.nan+np.nan*1j

    result = np.where((np.real(U2) > np.real(histo) + threshold * np.real(histostd)) |
                      (np.imag(U2) > np.imag(histo) + threshold * np.imag(histostd)) |
                      (np.real(U2) < np.real(histo) - threshold * np.real(histostd)) |
                      (np.imag(U2) < np.imag(histo) - threshold * np.imag(histostd)))

    for jj in range(len(result[0])):
        nu[result[0][jj], result[1][jj]] = np.nan
        nv[result[0][jj], result[1][jj]] = np.nan

    hu = nu[int(np.floor(m / 2)): -int(np.floor(m / 2)), int(np.floor(m / 2)): - int(np.floor(m / 2))]
    hv = nv[int(np.floor(m / 2)): -int(np.floor(m / 2)), int(np.floor(m / 2)): - int(np.floor(m / 2))]
    return hu, hv
