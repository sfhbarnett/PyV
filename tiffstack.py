import tifffile
import numpy as np
from skimage.transform import AffineTransform, warp

class tiffstack():

    """
    Container class for a stack of tiff images
    to do: handle both multipage tiffs as well as image sequences
    """

    def __init__(self, pathname=None):
        self.ims = None
        self.nfiles = 0
        self.minimum = 0
        self.maximum = np.inf
        self.width = 0
        self.height = 0
        self.dtype = None
        self.pathname = pathname
        if pathname is not None:
            self.load_info(pathname)
        self.transforms = []

    def load_info(self, pathname):
        self.ims = tifffile.TiffFile(pathname)
        self.width = self.ims.pages[0].shape[0]
        self.height = self.ims.pages[0].shape[1]
        self.dtype = self.ims.pages[0].dtype
        self.nfiles = len(self.ims.pages)

    def getimage(self, index):
        """
        Load in the image at index and store the min and max values
        :param index: which image in series to open
        :return: numpy array of image
        """
        image = self.ims.pages[index].asarray()
        self.minimum = image.min()
        self.maximum = image.max()
        return image

    def settransforms(self, xshift, yshift):
        """
        Once drift has been estimated, generate an affine transform for each image
        """
        self.transforms.append(AffineTransform(translation=[0, 0]))
        for i in range(self.nfiles-1):
            self.transforms.append(AffineTransform(translation=[xshift[i], yshift[i]]))

    def printtransforms(self):
        print(self.nfiles, len(self.transforms))
        for i in range(len(self.transforms)):
            print(self.transforms[i])

    def savedriftcorrected(self):
        outname = self.pathname[:-4] + 'DC.tif'
        with tifffile.TiffWriter(outname) as tif:
            for index in range(self.nfiles):
                image = self.getimage(index)
                transform = self.transforms[index]
                shifted = warp(image, transform, preserve_range=True)
                tif.save(np.int16(shifted))
        print('saved data')