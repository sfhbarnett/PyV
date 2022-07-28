import os

from PyV.PIV import PIV
from tiffstack import tiffstack
import matplotlib.pyplot as plt
from piv_filters import localfilt
from naninterp import naninterp
import matplotlib.cm as cm
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtGui import QIcon, QIntValidator, QDoubleValidator
import sys
import numpy as np
from main_gui import Ui_MainWindow
from superqt import QLabeledRangeSlider
import time
from skimage.transform import resize
import re

# pyuic6 - o main_gui.py - x main_gui.ui

class externalPIV(QtCore.QThread):
    progressstatus = QtCore.pyqtSignal(int)
    def __init__(self, x, y, u, v, windowsize, overlap, imstack):
        super(QtCore.QThread,self).__init__()
        self.u = u
        self.v = v
        self.x = x
        self.y = y
        self.windowsize = windowsize
        self.overlap = overlap
        self.imstack = imstack

    def run(self):
        for frame in range(self.imstack.nfiles-1):
            x, y, u, v = PIV(self.imstack.getimage(frame), self.imstack.getimage(frame+1), self.windowsize, self.overlap)
            u, v = localfilt(x, y, u, v, 2)
            u = naninterp(u)
            v = naninterp(v)
            self.x[:, :, frame] = x
            self.y[:, :, frame] = y
            self.u[:, :, frame] = u
            self.v[:, :, frame] = v
            self.progressstatus.emit(frame)
        return self.x,self.y,self.u,self.v


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.seticons()
        self.mincontrast = 0
        self.maxcontrast = 64000
        self.quiver = None
        self.mpl_toolbar = NavigationToolbar2QT(self.mplwidget.canvas, None)
        self.connectsignalsslots()
        self.stackslider.setRange(0, 10)
        self.currentimage = 0
        self.contrastslider.setRange(0, 200)
        self.contrastslider.setHandleLabelPosition(QLabeledRangeSlider.LabelPosition.LabelsBelow)
        self.contrastslider.setEdgeLabelMode(self.contrastslider.EdgeLabelMode.NoLabel)
        self.contrastslider.label_shift_x = 10
        self.maxdispspeed = int(self.maxspeedinput.text())
        self.maxspeedinput.setValidator(QDoubleValidator())
        self.windowsize = int(self.windowsizeinput.text())
        self.windowsizeinput.setValidator(QIntValidator())
        self.pixelsize = int(self.pixelsizeinput.text())
        self.pixelsizeinput.setValidator(QDoubleValidator())
        self.arrowscale = float(self.arrowscaleinput.text())
        self.arrowscaleinput.setValidator(QDoubleValidator())
        self.imagehandle = None
        self.progressbar = QtWidgets.QProgressBar()
        self.statusbar.addPermanentWidget(self.progressbar)

    def seticons(self):
        icon = QIcon("PyV/icons/pan.png")
        self.actionPan.setIcon(icon)
        icon = QIcon("PyV/icons/home.png")
        self.actionhome.setIcon(icon)
        icon = QIcon("PyV/icons/open.png")
        self.actionopen.setIcon(icon)
        icon = QIcon("PyV/icons/save.png")
        self.actionSave.setIcon(icon)
        icon = QIcon("PyV/icons/zoom2.png")
        self.actionzoom.setIcon(icon)
        icon = QIcon("PyV/icons/export.png")
        self.actionexport.setIcon(icon)
        icon = QIcon("PyV/icons/import.png")
        self.actionimport.setIcon(icon)

    def connectsignalsslots(self):
        self.actionzoom.triggered.connect(self.mpl_toolbar.zoom)
        self.actionPan.triggered.connect(self.mpl_toolbar.pan)
        self.actionhome.triggered.connect(self.mpl_toolbar.home)
        self.actionopen.triggered.connect(self.get_file)
        self.actionexport.triggered.connect(self.exportfields)
        self.actionimport.triggered.connect(self.importfields)
        self.runPIVbutton.clicked.connect(self.runPIV)
        self.alignmentbutton.clicked.connect(self.alignment)
        self.orientationbutton.clicked.connect(self.orientation)
        self.windowsizeinput.editingFinished.connect(self.setwindowsize)
        self.pixelsizeinput.editingFinished.connect(self.setPixelSize)
        self.maxspeedinput.editingFinished.connect(self.updatedisplayspeed)
        self.contrastslider.valueChanged.connect(self.update_contrast)
        self.stackslider.valueChanged.connect(self.move_through_stack)
        self.arrowscaleinput.editingFinished.connect(self.setarrowscale)
        self.actionSave.triggered.connect(self.savequiver)

    def get_file(self):
        self.mplwidget.canvas.axes.cla()
        self.filename = None
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', directory='~/Documents')
        if self.filename[0] != "":
            self.filename = self.filename[0]
            self.imstack = tiffstack(self.filename)
            self.imagehandle = self.mplwidget.canvas.axes.imshow(self.imstack.getimage(0))
            self.stackslider.setRange(0, self.imstack.nfiles - 1)
            self.stackslider.setValue(0)
            self.imagehandle.set_cmap('gray')
            self.mplwidget.canvas.fig.canvas.draw()
            self.contrastslider.setRange(0, self.imstack.maximum * 1.5)
            self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))
            self.stackpos.setText(str(self.currentimage + 1))

    def move_through_stack(self, value):
        """ Updates the current image in the viewport"""
        print(value)
        self.currentimage = value
        if self.imagehandle is not None:
            self.imagehandle.set_cmap('gray')
            self.stackpos.setText(str(self.currentimage+1))
            self.imagehandle.set_data(self.imstack.getimage(value))
            self.mplwidget.canvas.figure.canvas.draw_idle()
            # Check and update contrast
            if not self.autocontrast.isChecked():
                self.imagehandle.set_clim([self.mincontrast, self.maxcontrast])
            else:
                self.contrastslider.setRange(0, self.imstack.maximum*1.5)
                self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))
        if self.quiver is not None:
            self.makeQuiver()


    def update_contrast(self, value):
        self.mincontrast = value[0]
        self.maxcontrast = value[1]
        if self.imagehandle:
            self.imagehandle.set_clim(self.mincontrast, self.maxcontrast)
            self.mplwidget.canvas.fig.canvas.draw()

    def runPIV(self):
        overlap = 0.5
        self.x = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.y = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.u = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.v = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        # self.piv = externalPIV(self.x, self.y, self.u, self.v, self.windowsize, overlap, self.imstack)
        # self.piv.progressstatus.connect(self.updateprogress)
        # self.piv.start()
        start = time.time()
        for frame in range(self.imstack.nfiles-1):
            x, y, u, v = PIV(self.imstack.getimage(frame), self.imstack.getimage(frame+1), self.windowsize, overlap)
            u, v = localfilt(x, y, u, v, 2)
            u = naninterp(u)
            v = naninterp(v)
            self.x[:, :, frame] = x * self.pixelsize
            self.y[:, :, frame] = y * self.pixelsize
            self.u[:, :, frame] = u * self.pixelsize
            self.v[:, :, frame] = v * self.pixelsize
        print(time.time()-start)
        self.makeQuiver()


    def makeQuiver(self):
        if self.quiver is not None:
            self.quiver.remove()
        if self.currentimage != 0:
            frame = self.currentimage-1
            M = np.sqrt(self.u[:, :, frame]*self.u[:, :, frame]+self.v[:, :, frame]*self.v[:, :, frame])
            clipM = np.clip(M, 0, self.maxdispspeed)
            self.quiver = self.mplwidget.canvas.axes.quiver(self.x[:, :, frame]+1,
                                                            self.y[:, :, frame]+1,
                                                            self.u[:, :, frame]*self.arrowscale,
                                                            self.v[:, :, frame]*self.arrowscale, clipM,
                                                            scale_units='xy',scale=1, cmap=plt.cm.jet)
        else:
            self.quiver = self.mplwidget.canvas.axes.quiver([], [], [], [])
        self.mplwidget.canvas.fig.canvas.draw()

    def updatedisplayspeed(self):
        self.maxdispspeed = float(self.maxspeedinput.text())
        self.makeQuiver()

    def setPixelSize(self):
        value = float(self.pixelsizeinput.text())
        self.u /= self.pixelsize*self.windowsize
        self.v /= self.pixelsize*self.windowsize
        self.pixelsize = value
        self.u *= self.pixelsize*self.windowsize
        self.v *= self.pixelsize*self.windowsize
        self.makeQuiver()

    def setwindowsize(self):
        self.windowsize = int(self.windowsizeinput.text())

    def setarrowscale(self):
        self.arrowscale = float(self.arrowscaleinput.text())
        self.makeQuiver()

    def exportfields(self):
        path = '/Users/sbarnett/PycharmProjects/PyV/exportloc'
        fmt = '%d', '%d', '%1.3f', '%1.3f'
        for frame in range(self.x.shape[2]):
            array = np.vstack((self.x[:, :, frame].flatten('F'), self.y[:, :, frame].flatten('F'),
                               self.u[:, :, frame].flatten('F'), self.v[:, :, frame].flatten('F'))).T
            np.savetxt(path+'/'+str(frame)+'.txt', array, delimiter=',', fmt=fmt)
        self.statusbar.showMessage("PIV data exported to: "+path, 2000)

    def updateprogress(self, value):
        self.progressbar.setValue(100 / self.imstack.nfiles * value)

    def savequiver(self):
        extent = self.mplwidget.canvas.axes.get_window_extent().transformed(self.mplwidget.canvas.fig.dpi_scale_trans.inverted())
        for frame in range(1, self.imstack.nfiles):
            self.move_through_stack(frame)
            self.mplwidget.canvas.fig.savefig('/Users/sbarnett/PycharmProjects/PyV/exportloc/'+str(frame)+'.png',
                                              bbox_inches=extent)

    def alignment(self):
        frame = self.currentimage
        averageU = np.mean(self.u[:,:,frame])
        averageV = np.mean(self.v[:,:,frame])
        meanvector = np.array([averageU,averageV])
        upper = np.vstack((self.u[:,:,frame].flatten('F'),self.v[:,:,frame].flatten('F'))).T@meanvector
        lower = np.sqrt(self.u[:, :, frame].flatten('F')**2+self.v[:, :, frame].flatten('F')**2) * np.sqrt(meanvector[0]**2+meanvector[1]**2)
        value = np.divide(upper,lower)
        map = np.reshape(value,(63,63))
        map = resize(map,(self.u.shape[0]*self.windowsize//2,self.u.shape[1]*self.windowsize//2),order=0,preserve_range=True)
        self.imagehandle.set_data(map)
        self.imagehandle.set_cmap('viridis')
        self.imagehandle.set_clim(-1, 1)
        self.mplwidget.canvas.fig.canvas.draw()


    def orientation(self):
        rho = np.sqrt(self.u[:, :, self.currentimage]**2+self.v[:, :, self.currentimage]**2)
        phi = np.arctan2(self.v[:, :, self.currentimage], self.u[:, :, self.currentimage])
        theta = resize(np.rad2deg(phi), (self.u.shape[0]*self.windowsize//2,self.u.shape[1]*self.windowsize//2), order=1, preserve_range=True)
        self.imagehandle.set_data(theta+180)
        self.imagehandle.set_clim(0, 360)
        self.imagehandle.set_cmap('hsv')
        self.mplwidget.canvas.fig.canvas.draw()

    def importfields(self):
        """ import files that represents fields (Mx4) where M is the number of vectors. Col 0 is x, Col 1 is y,
        Col 2 is u, Col3 is v"""
        path = '/Users/sbarnett/PycharmProjects/PyV/exportloc/'
        files = os.listdir(path)
        files.remove('.DS_Store')
        counter = 0
        self.x = np.zeros((63,63,9))
        self.y = np.zeros((63,63,9))
        self.u = np.zeros((63,63,9))
        self.v = np.zeros((63,63,9))
        files = self.natsort(files)
        for file in files:
            if file[-4:] == '.txt':
                array = np.loadtxt(os.path.join(path,file), delimiter=',')
                self.x[:, :, counter] = np.reshape(array[:,0],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.y[:, :, counter] = np.reshape(array[:,1],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.u[:, :, counter] = np.reshape(array[:,2],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.v[:, :, counter] = np.reshape(array[:,3],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                counter += 1
        self.makeQuiver()

    def natsort(self,tosort):
        """Natural sort file strings"""
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(tosort, key=alphanum_key)



def main():
    #Initialise GUI
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    app.setQuitOnLastWindowClosed(True)

    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()