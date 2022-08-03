import numpy as np
from numba import njit
import pyfftw
import multiprocessing

@njit(nopython=True)
def intpeak(x1, y1, R, Rxm1, Rxp1, Rym1, Ryp1, N):
    M = N
    x01 = x1 + ((np.log(Rxm1) - np.log(Rxp1)) / ((2 * np.log(Rxm1)) - (4 * np.log(R)) + (2 * np.log(Rxp1))))
    y01 = y1 + ((np.log(Rym1) - np.log(Ryp1)) / ((2 * np.log(Rym1)) - (4 * np.log(R)) + (2 * np.log(Ryp1))))
    x0 = x01-(M-1)
    y0 = y01-(N-1)
    return x0, y0


def nextpow2(x):
    n = 0
    while 2**n < x:
        n += 1
    return n


def weight(siz, dev):
    w = 1-np.power(np.cos(np.pi*np.array([x for x in range(siz)])/(siz-1)), dev)
    w2 = 1-np.power(np.cos(np.pi*np.array([x for x in range(siz)])/(siz-1)), dev)
    w = np.outer(w.T, w2)
    return w


def xcorrf2(crop1, crop2, mf, nf):
    #pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
    ma = crop1.shape[0]
    na = crop1.shape[1]
    mb = crop2.shape[0]
    nb = crop2.shape[1]

    b = np.conj(crop2[mb::-1, nb::-1])
    #at = np.fft.fft2(b, s=(mf, nf))
    at = pyfftw.builders.fft2(b, s=(mf, nf),threads=8)
    #bt = np.fft.fft2(crop1, s=(mf, nf))
    bt = pyfftw.builders.fft2(crop1, s=(mf, nf),threads=8)
    #c = np.fft.ifft2(np.multiply(at, bt))
    test = np.multiply(at,bt)
    c = pyfftw.builders.ifft2(test,threads=8)
    c = np.real(c)
    output = c[:ma + mb-1, :na + nb-1]
    return output
