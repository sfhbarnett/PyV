from PyV.PIV import PIV
from tiffstack import TiffStack
import matplotlib.pyplot as plt
import numpy as np
from piv_filters import localfilt


path = '/Users/sbarnett/Documents/PIVData/fatima/ForSam/C1-20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi - 20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi #21.tif'

tf = TiffStack(path)
image1 = tf.getimage(0)
image2 = tf.getimage(1)

# plt.imshow(tf.getimage(0))
# plt.show()
overlap = 0.5
windowsize = 16

x, y, u, v = PIV(image1, image2, windowsize, overlap)

u, v = localfilt(x, y, u, v, 2)

plt.quiver(x, y, u, v, scale_units='dots', width=0.001, headlength=3, headwidth=2)
plt.axis('equal')
plt.show()
