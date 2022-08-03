import numpy as np
from utils import intpeak, weight, xcorrf2, nextpow2



def PIV(image1,image2,windowsize,overlap):
    stdim1 = np.std(image1)
    IN = np.zeros(image1.shape)
    ci1 = 0
    cj1 = 0
    sx = len([x for x in range(0, image1.shape[0]-windowsize+1, int((1-overlap)*windowsize))])
    sy = len([y for y in range(0, image1.shape[1]-windowsize+1, int((1-overlap)*windowsize))])
    mf = 2 ** nextpow2(windowsize + windowsize)
    nf = mf
    BiCor = np.divide(xcorrf2(np.ones((windowsize, windowsize)), np.ones((windowsize, windowsize)),mf,nf), windowsize**2)
    BiCor = BiCor/np.max(BiCor)
    x = np.zeros((sx, sy))
    y = np.zeros((sx, sy))
    u = np.zeros((sx, sy))
    v = np.zeros((sx, sy))
    dt = 1
    w = weight(windowsize, 20)
    windowsizesquare = windowsize**2

    for jj in range(0, image1.shape[0]-windowsize//2, int((1-overlap)*windowsize)):
        for ii in range(0, image1.shape[1]-windowsize//2, int((1-overlap)*windowsize)):
            if IN[int(jj+windowsize/2), int(ii+windowsize/2)] != 1:
                crop1 = image1[jj:jj+windowsize, ii:ii+windowsize]
                crop2 = image2[jj:jj+windowsize, ii:ii+windowsize]
                stdcrop1 = np.std(crop1)
                stdcrop2 = np.std(crop2)
                if stdcrop1 < 0.2*stdim1:
                    stdcrop1 = np.nan
                if stdcrop2 < 0.2*stdim1:
                    stdcrop2 = np.nan

                crop1 = crop1-np.mean(crop1)
                crop2 = crop2-np.mean(crop2)
                crop1 = np.multiply(crop1, w)
                crop2 = np.multiply(crop2, w)

                if np.isnan(stdcrop1) != 1 and np.isnan(stdcrop2) != 1:
                    output = xcorrf2(crop1, crop2,mf,nf)
                    output = output/(windowsizesquare*stdcrop1*stdcrop2)
                    output = np.divide(output, BiCor)
                    maxima = np.where(output == np.nanmax(output[int(0.5*windowsize+1):int(1.5*windowsize-3), int(0.5*windowsize+1):int(1.5*windowsize-3)]))
                    if len(maxima) > 1:
                        maxima = list(zip(maxima[0], maxima[1]))
                        maxima = maxima[0]
                    y1 = maxima[0]
                    x1 = maxima[1]
                    x0, y0 = intpeak(x1, y1, output[y1, x1], output[y1, x1-1], output[y1, x1+1], output[y1-1, x1], output[y1+1, x1], windowsize)
                    x[cj1, ci1] = int(windowsize/2+ii-1)
                    y[cj1, ci1] = int(windowsize / 2 + jj - 1)
                    u[cj1, ci1] = -x0/dt
                    v[cj1, ci1] = y0/dt
            ci1 += 1
        cj1 += 1
        ci1 = 0

    return x, y, u, v

def PIVvect(image1,image2,windowsize,overlap):
    stdim1 = np.std(image1)
    IN = np.zeros(image1.shape)
    ci1 = 0
    cj1 = 0
    sx = len([x for x in range(0, image1.shape[0]-windowsize+1, int((1-overlap)*windowsize))])
    sy = len([y for y in range(0, image1.shape[1]-windowsize+1, int((1-overlap)*windowsize))])
    mf = 2 ** nextpow2(windowsize + windowsize)
    nf = mf
    BiCor = np.divide(xcorrf2(np.ones((windowsize, windowsize)), np.ones((windowsize, windowsize)),mf,nf), windowsize**2)
    BiCor = BiCor/np.max(BiCor)
    x = np.zeros((sx, sy))
    y = np.zeros((sx, sy))
    u = np.zeros((sx, sy))
    v = np.zeros((sx, sy))
    dt = 1
    w = weight(windowsize, 20)
    windowsizesquare = windowsize**2
    im1 = np.lib.stride_tricks.sliding_window_view(image1[:,:],(windowsize,windowsize))[::int((1-overlap)*windowsize),::int((1-overlap)*windowsize)]
    im2 = np.lib.stride_tricks.sliding_window_view(image2[:,:],(windowsize,windowsize))[::int((1-overlap)*windowsize),::int((1-overlap)*windowsize)]
    im1 = im1.reshape(-1, *im1.shape[2:])
    im2 = im2.reshape(-1, *im2.shape[2:])
    std1 = im1.std(axis=(1,2))
    std2 = im2.std(axis=(1,2))
    std1[np.where(std1 < stdim1*0.2)] = np.nan
    std2[np.where(std2 < stdim1 * 0.2)] = np.nan
    im1 = im1 - im1.mean(axis=(1, 2), keepdims=True)
    im2 = im2 - im2.mean(axis=(1, 2), keepdims=True)
    im1 = np.multiply(im1, w)
    im2 = np.multiply(im2, w)

    #perform cross crorrelation
    t = time.time()
    at = np.fft.fft2(np.conj(im2[:, windowsize::-1, windowsize::-1]), s=(mf, nf))
    c = np.fft.ifft2(np.multiply(at, np.fft.fft2(im1, s=(mf, nf))))
    c = np.real(c)
    print(time.time()-t)
    cc = c[:, :windowsize + windowsize - 1, :windowsize + windowsize - 1]
    cc = cc / (windowsizesquare * std1[:, np.newaxis, np.newaxis] * std2[:, np.newaxis, np.newaxis])
    cc = cc / BiCor[np.newaxis, :]
    maxima = np.array([np.where(slice == np.nanmax(slice[int(0.5 * windowsize + 1):int(1.5 * windowsize - 3),
                                          int(0.5 * windowsize + 1):int(1.5 * windowsize - 3)])) for slice in cc])
    #might be error where multiple maxima found, have to get first element from each entry
    y1 = maxima[:, 0]
    x1 = maxima[:, 1]
    locs = np.array([intpeak(x2, y2, output[y2, x2], output[y2, x2 - 1], output[y2, x2 + 1], output[y2 - 1, x2],
                             output[y2 + 1, x2], windowsize) for x2,y2,output in zip(x1,y1,cc)]).squeeze()

    u = np.reshape(locs[:,0],[63,63])/dt
    v = np.reshape(locs[:,1],[63,63])/dt
    x = np.tile(np.arange(0,63),(63,1))
    y = x.T
    return x,y,u,v



if __name__ == '__main__':
    import cProfile
    from tiffstack import tiffstack
    import matplotlib.pyplot as plt

    path = '/Users/sbarnett/Documents/PIVData/fatima/ForSam/monolayer2/' \
           'C1-20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi -' \
           ' 20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi #21.tif'
    tf = tiffstack(path)
    #cProfile.run('PIVvect(tf.getimage(0),tf.getimage(1),32,0.5)')
    import time
    start = time.time()
    PIV(tf.getimage(0), tf.getimage(1), 32, 0.5)
    print(time.time()-start)