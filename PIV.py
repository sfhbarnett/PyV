import numpy as np
from utils import intpeak, weight, xcorrf2


def PIV(image1,image2,windowsize,overlap):
    stdim1 = np.std(image1)
    IN = np.zeros(image1.shape)
    ci1 = 0
    cj1 = 0
    sx = len([x for x in range(0, image1.shape[0]-windowsize+1, int((1-overlap)*windowsize))])
    sy = len([y for y in range(0, image1.shape[1]-windowsize+1, int((1-overlap)*windowsize))])
    BiCor = np.divide(xcorrf2(np.ones((windowsize, windowsize)), np.ones((windowsize, windowsize))), windowsize**2)
    BiCor = BiCor/np.max(BiCor)
    x = np.zeros((sx, sy))
    y = np.zeros((sx, sy))
    u = np.zeros((sx, sy))
    v = np.zeros((sx, sy))
    dt = 1
    w = weight(windowsize, 20)

    for jj in range(0, image1.shape[0]-windowsize+1, int((1-overlap)*windowsize)):
        for ii in range(0, image1.shape[1]-windowsize+1, int((1-overlap)*windowsize)):
            if cj1 == 65 and ci1 == 69:
                a = 'break'
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
                    output = xcorrf2(crop1, crop2)
                    output = output/(windowsize**2*stdcrop1*stdcrop2)
                    output = np.divide(output, BiCor)
                    maxima = np.where(output == np.nanmax(output[int(0.5*windowsize+1):int(1.5*windowsize-3), int(0.5*windowsize+1):int(1.5*windowsize-3)]))
                    if len(maxima) > 1:
                        maxima = list(zip(maxima[0], maxima[1]))
                        maxima = maxima[0]
                    y1 = maxima[0]
                    x1 = maxima[1]
                    x0,y0 = intpeak(x1, y1, output[y1, x1], output[y1, x1-1], output[y1, x1+1], output[y1-1, x1], output[y1+1, x1], windowsize)
                    # R2 = output
                    # R2[y1-3:y1+3,x1-3:x1+3] = np.nan
                    # p2_y2,p2_x2 = np.where(R2 == np.amax(R2))
                    x[cj1, ci1] = int(windowsize/2+ii-1)
                    y[cj1, ci1] = int(windowsize / 2 + jj - 1)
                    u[cj1, ci1] = -x0/dt
                    v[cj1, ci1] = y0/dt
            ci1 += 1
        cj1 += 1
        ci1 = 0

    return x, y, u, v