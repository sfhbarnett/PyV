import os
from PyV.PIV import PIV
from tiffstack import tiffstack
from piv_filters import localfilt
from naninterp import naninterp
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6 import QtCore, QtWidgets
import matplotlib.pyplot as plt
from PyQt6.QtGui import QIcon, QIntValidator, QDoubleValidator
from main_gui import Ui_MainWindow
from superqt import QLabeledDoubleRangeSlider
import multiprocessing
from skimage.transform import resize
import re
import h5py
import sys
import numpy as np

# pyuic6 - o main_gui.py - x main_gui.ui

class WorkerSignals(QtCore.QObject):

    finished = QtCore.pyqtSignal()
    error = QtCore.pyqtSignal(tuple)
    result = QtCore.pyqtSignal(tuple)
    progress = QtCore.pyqtSignal(int)


class Worker(QtCore.QRunnable):

    def __init__(self, *args, **kwargs):
        super(Worker, self).__init__()
        self.imstack = args[0]
        self.windowsize = args[1]
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.progress = 0

    @QtCore.pyqtSlot()
    def run(self):
        with multiprocessing.Pool() as pool:
            for frame in range(self.imstack.nfiles-1):
                pool.apply_async(pivwrapper, args=(self.imstack.getimage(frame), self.imstack.getimage(frame+1), self.windowsize, frame),callback=self.emitresults)
            pool.close()
            pool.join()

    @QtCore.pyqtSlot()
    def emitresults(self, results):
        self.progress+=1
        self.signals.result.emit(results)
        self.signals.progress.emit(self.progress)

def pivwrapper(image1, image2, windowsize, frame):
    print("working frame " + str(frame))
    overlap = 0.5
    x, y, u, v = PIV(image1, image2, windowsize, overlap)
    u, v = localfilt(x, y, u, v, 2)
    u = naninterp(u)
    v = naninterp(v)
    return x, y, u, v, frame


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.setWindowTitle("PyV - v0.1")
        self.seticons()
        self.piv = Pivdata()
        self.mincontrast = 0
        self.maxcontrast = 64000
        self.quiver = None
        self.imagehandle = plt.imshow(np.array([[]]))
        self.imagehandle.set_cmap('gray')
        self.extent = [0, 1000, 0, 1000]
        self.imstack = None
        self.showQuiver = False
        self.mpl_toolbar = NavigationToolbar2QT(self.mplwidget.canvas, None)
        self.connectsignalsslots()
        self.stackslider.setRange(0, 10)
        self.currentimage = 0
        self.contrastslider.setRange(0, 200)
        self.contrastslider.setHandleLabelPosition(QLabeledDoubleRangeSlider.LabelPosition.LabelsBelow)
        self.contrastslider.setEdgeLabelMode(self.contrastslider.EdgeLabelMode.NoLabel)
        self.contrastslider.label_shift_x = 10
        self.contrastslider.setSingleStep(1)
        self.contrastslider.setDecimals(0)
        self.maxdispspeed = int(self.maxspeedinput.text())
        self.maxspeedinput.setValidator(QDoubleValidator())
        self.windowsize = int(self.windowsizeinput.text())
        self.windowsizeinput.setValidator(QIntValidator())
        self.pixelsize = int(self.pixelsizeinput.text())
        self.pixelsizeinput.setValidator(QDoubleValidator())
        self.arrowscale = float(self.arrowscaleinput.text())
        self.arrowscaleinput.setValidator(QDoubleValidator())
        self.timeinterval = float(self.timeintervalinput.text())
        self.timeintervalinput.setValidator(QDoubleValidator())
        self.centerXinput.setValidator(QDoubleValidator())
        self.centerYinput.setValidator(QDoubleValidator())
        self.tableview = TableView(self.table)
        self.colormaps = ['jet', 'viridis', 'plasma', 'inferno']
        self.solids = ['black', 'green', 'red', 'cyan', 'magenta', 'yellow']
        self.solidcode = {'black': 'k', 'green': 'g', 'red': 'r', 'cyan': 'c', 'magenta': 'm', 'yellow': 'y'}
        self.fieldcolormapinput.addItems(self.colormaps+self.solids)
        self.quivercmap = 'jet'
        self.quivercolor = None
        self.progressbar = QtWidgets.QProgressBar()
        self.statusbar.addPermanentWidget(self.progressbar)
        self.alignmentbutton.setStyleSheet("")
        self.threadpool = QtCore.QThreadPool()
        self.setFocus()

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
        # Input output
        self.actionopen.triggered.connect(self.get_file)
        self.actionexport.triggered.connect(self.exportfields)
        self.actionimport.triggered.connect(self.importfields)
        self.actionSave.triggered.connect(self.savequiver)
        self.actionExport_as_hdf5.triggered.connect(self.exporth5py)
        self.actionImport_fields_from_hdf5.triggered.connect(self.importh5py)
        # PIV
        self.runPIVbutton.clicked.connect(self.runPIV)
        self.windowsizeinput.editingFinished.connect(self.setwindowsize)
        self.pixelsizeinput.editingFinished.connect(self.setPixelSize)
        self.maxspeedinput.editingFinished.connect(self.updatedisplayspeed)
        self.arrowscaleinput.editingFinished.connect(self.setarrowscale)
        # Analysis
        self.alignmentbutton.clicked.connect(self.alignment)
        self.orientationbutton.clicked.connect(self.orientation)
        self.timeintervalinput.editingFinished.connect(self.settimeinterval)
        self.vrmsbutton.clicked.connect(self.calculatevrms)
        self.linearorderbutton.clicked.connect(self.calculateLinearOrderParameter)
        self.linearisebutton.clicked.connect(self.lineariseField)
        self.centerXinput.editingFinished.connect(self.setCenterX)
        self.centerYinput.editingFinished.connect(self.setCenterY)
        self.rotationalorderbutton.clicked.connect(self.rotationalorder)
        # Display
        self.actionzoom.triggered.connect(self.mpl_toolbar.zoom)
        self.actionPan.triggered.connect(self.mpl_toolbar.pan)
        self.actionhome.triggered.connect(self.mpl_toolbar.home)
        self.contrastslider.valueChanged.connect(self.update_contrast)
        self.stackslider.valueChanged.connect(self.move_through_stack)
        self.fieldcolormapinput.currentIndexChanged.connect(self.changecmap)
        self.showimagecheck.stateChanged.connect(self.showImage)
        self.showlinearfield.stateChanged.connect(self.showlinear)

    def get_file(self):
        """"
        Load in image stack of raw data
        """
        self.mplwidget.canvas.axes.cla()
        self.filename = None
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', directory='~/Documents')
        if self.filename[0] != "":
            self.filename = self.filename[0]
            self.imstack = tiffstack(self.filename)
            self.extent = [0, self.imstack.width, 0, self.imstack.height]
            self.imagehandle = self.mplwidget.canvas.axes.imshow(self.imstack.getimage(0), extent=self.extent)
            self.stackslider.setRange(0, self.imstack.nfiles - 1)
            self.stackslider.setValue(0)
            self.imagehandle.set_cmap('gray')
            self.mplwidget.canvas.fig.canvas.draw()
            self.contrastslider.setRange(0, self.imstack.maximum * 1.5)
            self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))
            self.stackpos.setText(str(self.currentimage + 1))

    def move_through_stack(self, value):
        """ Updates the current image in the viewport"""
        self.currentimage = value
        self.showImage()
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
        self.progress = 0
        self.piv.x = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.y = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.u = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.v = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        # Create worker thread to perform multiprocessing
        self.threadpool = QtCore.QThreadPool()
        worker = Worker(self.imstack, self.windowsize)
        worker.signals.result.connect(self.assignpivresult)
        worker.signals.progress.connect(self.updateprogress)
        self.threadpool.start(worker)
        self.showQuiver = True
        self.makeQuiver()

    def assignpivresult(self, result):
        x, y, u, v, frame = result
        self.piv.x[:, :, frame] = x * self.pixelsize
        self.piv.y[:, :, frame] = y * self.pixelsize
        self.piv.u[:, :, frame] = u * self.pixelsize / self.timeinterval
        self.piv.v[:, :, frame] = v * self.pixelsize / self.timeinterval
        self.progress+=1

    def makeQuiver(self):
        if self.quiver is not None:
            self.quiver.remove()
        if self.currentimage != 0 and self.showQuiver is True:
            frame = self.currentimage-1
            if self.showlinearfield.isChecked():
                u = self.piv.ru[:, :, frame]
                v = self.piv.rv[:, :, frame]
            else:
                u = self.piv.u[:, :, frame]
                v = self.piv.v[:, :, frame]
            M = np.sqrt(u**2 + v**2)
            clipM = np.clip(M, 0, self.maxdispspeed)
            if self.quivercmap is not None:
                self.quiver = self.mplwidget.canvas.axes.quiver(self.piv.x[:, :, frame]+1,
                                                                self.piv.y[:, :, frame]+1,
                                                                u*self.arrowscale,
                                                                v*self.arrowscale, clipM,
                                                                scale_units='xy', scale=1,
                                                                cmap=self.quivercmap)
            else:
                self.quiver = self.mplwidget.canvas.axes.quiver(self.piv.x[:, :, frame]+1,
                                                                self.piv.y[:, :, frame]+1,
                                                                u*self.arrowscale,
                                                                v*self.arrowscale,
                                                                color=self.quivercolor,
                                                                scale_units='xy', scale=1)
        else:
            self.quiver = self.mplwidget.canvas.axes.quiver([], [], [], [])
        self.imagehandle.set_extent(self.extent)
        self.mplwidget.canvas.fig.canvas.draw()

    def updatedisplayspeed(self):
        self.maxdispspeed = float(self.maxspeedinput.text())
        self.makeQuiver()
        self.setFocus()

    def setPixelSize(self):
        value = float(self.pixelsizeinput.text())
        self.piv.u /= self.pixelsize*self.windowsize
        self.piv.v /= self.pixelsize*self.windowsize
        self.pixelsize = value
        self.piv.u *= self.pixelsize*self.windowsize
        self.piv.v *= self.pixelsize*self.windowsize
        self.makeQuiver()
        self.setFocus()

    def setCenterX(self):
        value = float(self.centerXinput.text())
        self.piv.centerx = value
        self.setFocus()

    def setCenterY(self):
        value = float(self.centerYinput.text())
        self.piv.centery = value
        self.setFocus()

    def setwindowsize(self):
        self.windowsize = int(self.windowsizeinput.text())
        self.setFocus()

    def setarrowscale(self):
        """
        Sets a scaling factor for quiver arrows, doesn't affect underlying data
        :return:
        """
        self.arrowscale = float(self.arrowscaleinput.text())
        self.makeQuiver()
        self.setFocus()

    def settimeinterval(self):
        value = float(self.timeintervalinput.text())
        self.piv.u *= self.timeinterval
        self.piv.v *= self.timeinterval
        self.timeinterval = value
        self.piv.u /= self.timeinterval
        self.piv.v /= self.timeinterval
        self.makeQuiver()
        self.setFocus()

    def changecmap(self):
        cmap = self.fieldcolormapinput.currentText()
        if cmap in self.colormaps:
            self.quivercmap = cmap
            self.quivercolor = None
        else:
            self.quivercmap = None
            self.quivercolor = cmap
        self.makeQuiver()

    def exportfields(self):
        """
        Export to csv
        :return:
        """
        path = '/Users/sbarnett/PycharmProjects/PyV/exportloc'
        fmt = '%d', '%d', '%1.3f', '%1.3f'
        for frame in range(self.piv.x.shape[2]):
            array = np.vstack((self.piv.x[:, :, frame].flatten('F'), self.piv.y[:, :, frame].flatten('F'),
                               self.piv.u[:, :, frame].flatten('F'), self.piv.v[:, :, frame].flatten('F'))).T
            np.savetxt(path+'/'+str(frame)+'.txt', array, delimiter=',', fmt=fmt)
        self.statusbar.showMessage("PIV data exported to: "+path, 2000)

    def exporth5py(self):
        """
        Exports fields and relevent information to hdf5
        :return:
        """
        filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', directory='~/Documents')
        if filename[0] != "":
            hf = h5py.File(filename[0], 'w')
            data = hf.create_group('data')
            data.create_dataset('xfield', data=self.piv.x)
            data.create_dataset('yfield', data=self.piv.y)
            data.create_dataset('ufield', data=self.piv.u)
            data.create_dataset('vfield', data=self.piv.v)
            params = hf.create_group('params')
            params.create_dataset('windowsize', data=self.piv.windowsize)
            params.create_dataset('pixelsize', data=self.piv.pixelsize)
            params.create_dataset('timeinterval', data=self.piv.timeinterval)
            params.create_dataset('width', data=self.piv.width)
            params.create_dataset('height', data=self.piv.height)
            params.create_dataset('nfields', data=self.piv.nfields)
            hf.close()
            self.statusbar.showMessage("PIV data exported to: " + filename[0], 2000)

    def importh5py(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', directory='~/Documents')
        with h5py.File(filename[0], 'r') as hf:
            data = hf.get('data')
            self.piv.x = np.array(data.get('xfield'))
            self.piv.y = np.array(data.get('yfield'))
            self.piv.u = np.array(data.get('ufield'))
            self.piv.v = np.array(data.get('vfield'))
            params = hf.get('params')
            self.piv.windowsize = int(np.array(params.get('windowsize')))
            self.piv.pixelsize = float(np.array(params.get('pixelsize')))
            self.piv.timeinterval = float(np.array(params.get('timeinterval')))
            self.piv.width = int(np.array(params.get('width')))
            self.piv.height = int(np.array(params.get('height')))
            self.piv.nfields = int(np.array(params.get('nfields')))
            self.showQuiver = True
            self.statusbar.showMessage("PIV data imported from: " + filename[0], 2000)
            self.pixelsizeinput.setText(str(self.piv.pixelsize))
            self.windowsizeinput.setText(str(self.piv.windowsize))
            self.timeintervalinput.setText(str(self.piv.timeinterval))


    def updateprogress(self, value):
        self.progressbar.setValue((value / (self.imstack.nfiles-1))*100)

    def savequiver(self):
        """
        Saves the current image on the screen
        :return:
        """
        extent = self.mplwidget.canvas.axes.get_window_extent().transformed(self.mplwidget.canvas.fig.dpi_scale_trans.inverted())
        for frame in range(1, self.imstack.nfiles):
            self.move_through_stack(frame)
            self.mplwidget.canvas.fig.savefig('/Users/sbarnett/PycharmProjects/PyV/exportloc/'+str(frame)+'.png',
                                              bbox_inches=extent)

    def alignment(self):
        """
        Create an alignment map where the pixels are colored based on how well they align with the mean vector
        :return:
        """
        if self.currentimage != 0:
            frame = self.currentimage-1
            averageU = np.mean(self.piv.u[:, :, frame])
            averageV = np.mean(self.piv.v[:, :, frame])
            meanvector = np.array([averageU, averageV])
            upper = np.vstack((self.piv.u[:, :, frame].flatten('F'), self.piv.v[:, :, frame].flatten('F'))).T@meanvector
            lower = np.sqrt(self.piv.u[:, :, frame].flatten('F')**2 + self.piv.v[:, :, frame].flatten('F')**2)\
                    * np.sqrt(meanvector[0]**2+meanvector[1]**2)
            value = np.divide(upper, lower)
            alignmap = np.reshape(value, (63, 63))
            alignmap = resize(alignmap, (self.piv.u.shape[0]*self.windowsize//2,
                                         self.piv.u.shape[1]*self.windowsize//2), order=0, preserve_range=True)
            self.imagehandle.set_data(np.flipud(alignmap))
            self.contrastslider.setRange(-1, 1)
            self.contrastslider.setValue((-1, 1))
            self.contrastslider.setSingleStep(0.01)
            self.contrastslider.setDecimals(2)
            self.contrastslider.label_shift_x = 10
            self.extent = [self.windowsize/4, self.imstack.width-self.windowsize/4,
                           self.windowsize/4, self.imstack.width-self.windowsize/4]
            self.imagehandle.set_extent(self.extent)
            self.imagehandle.set_cmap('viridis')
            self.imagehandle.set_clim(-1, 1)
            self.mplwidget.canvas.fig.canvas.draw()

    def orientation(self):
        """
        Create an orientation map where the image is colored based on the direction of the vector
        :return: None
        """
        if self.currentimage != 0:
            frame = self.currentimage-1
            rho = np.sqrt(self.piv.u[:, :, frame]**2+self.piv.v[:, :, frame]**2)
            phi = np.arctan2(self.piv.v[:, :, frame], self.piv.u[:, :, frame])
            theta = resize(np.rad2deg(phi), (self.piv.u.shape[0]*self.windowsize//2,
                                             self.piv.u.shape[1]*self.windowsize//2),
                           order=0, preserve_range=True)
            self.imagehandle.set_data(np.flipud(theta.T)+180)
            self.imagehandle.set_clim(0, 360)
            self.extent = [self.windowsize/4, self.imstack.width-self.windowsize/4,
                           self.windowsize/4, self.imstack.width-self.windowsize/4]
            self.imagehandle.set_extent(self.extent)
            self.imagehandle.set_cmap('hsv')
            self.mplwidget.canvas.fig.canvas.draw()

    def importfields(self):
        """ import files that represents fields (Mx4) where M is the number of vectors. Col 0 is x, Col 1 is y,
        Col 2 is u, Col3 is v"""
        path = '/Users/sbarnett/PycharmProjects/PyV/exportloc/'
        files = os.listdir(path)
        files.remove('.DS_Store')
        counter = 0
        self.piv.x = np.zeros((63, 63, 9))
        self.piv.y = np.zeros((63, 63, 9))
        self.piv.u = np.zeros((63, 63, 9))
        self.piv.v = np.zeros((63, 63, 9))
        files = self.natsort(files)
        for file in files:
            if file[-4:] == '.txt':
                array = np.loadtxt(os.path.join(path, file), delimiter=',')
                self.piv.x[:, :, counter] = np.reshape(array[:, 0], [int(np.max(array[:, 0]/(array[0, 0]+1)))+1,
                                                                     int(np.max(array[:, 1]/(array[0, 1]+1)))+1])
                self.piv.y[:, :, counter] = np.reshape(array[:, 1], [int(np.max(array[:, 0]/(array[0, 0]+1)))+1,
                                                                     int(np.max(array[:, 1]/(array[0, 1]+1)))+1])
                self.piv.u[:, :, counter] = np.reshape(array[:, 2], [int(np.max(array[:, 0]/(array[0, 0]+1)))+1,
                                                                     int(np.max(array[:, 1]/(array[0, 1]+1)))+1])
                self.piv.v[:, :, counter] = np.reshape(array[:, 3], [int(np.max(array[:, 0]/(array[0, 0]+1)))+1,
                                                                     int(np.max(array[:, 1]/(array[0, 1]+1)))+1])
                counter += 1
        self.piv.width = int(np.max(array[:, 0] / (array[0, 0] + 1))) + 1
        self.piv.height = int(np.max(array[:, 1] / (array[0, 1] + 1))) + 1
        self.piv.nfields = counter
        self.showQuiver = True
        self.makeQuiver()
        self.statusbar.showMessage("PIV data imported from: "+path, 2000)

    def calculatevrms(self):
        """Calculate root mean square velocity of the vectorfield"""
        self.vrms = []
        for frame in range(self.piv.u.shape[2]):
            velocities = np.sqrt(self.piv.u[:, :, frame]**2 + self.piv.v[:, :, frame]**2)
            self.vrms.append(np.sqrt(np.nanmean(velocities**2)))
        self.tableview.addData(self.vrms, 'Vrms')

    def rotationalorder(self):
        pass

    def lineariseField(self):
        self.piv.lineariseField()
        self.makeQuiver()

    def showImage(self):
        if self.imstack is not None and self.showimagecheck.isChecked():
            self.stackpos.setText(str(self.currentimage+1))
            self.imagehandle.set_data(self.imstack.getimage(self.currentimage))
            self.imagehandle.set_extent(self.extent)
            # Check and update contrast
            if not self.autocontrast.isChecked():
                self.imagehandle.set_clim([self.mincontrast, self.maxcontrast])
            else:
                self.contrastslider.setRange(0, self.imstack.maximum*1.5)
                self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))
        else:
            self.imagehandle.set_data(np.array([[], []]))
        self.imagehandle.set_cmap('gray')
        self.mplwidget.canvas.figure.canvas.draw_idle()

    def showlinear(self):
        self.makeQuiver()

    def calculateLinearOrderParameter(self):
        """
        Linear order parameter describes how well aligned all the vectors are in a particular direction.
        A value of 1 indicates perfect alignment and a value of 0 indicates chaos
        """
        self.LOP = []
        for frame in range(self.piv.u.shape[2]):
            mean_u = np.nanmean(self.piv.u[:, :, frame])
            mean_v = np.nanmean(self.piv.v[:, :, frame])
            mean_uv2 = np.nanmean(self.piv.u[:, :, frame]**2 + self.piv.v[:, :, frame]**2)
            self.LOP.append((mean_u**2+mean_v**2)/mean_uv2)
        self.tableview.addData(self.LOP, 'LinearOrder')

    def natsort(self, tosort):
        """Natural sort file strings"""
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(tosort, key=alphanum_key)

    def keyPressEvent(self, event):
        # Left and right for moving through stack
        if self.imagehandle is not None or self.quiver is not None:
            if event.key() == QtCore.Qt.Key_S or event.key() == QtCore.Qt.Key_Right:
                if 0 <= self.stackslider.value() < self.imstack.nfiles-1:
                    self.stackslider.setValue(self.stackslider.value()+1)
            elif event.key() == QtCore.Qt.Key_A or event.key() == QtCore.Qt.Key_Left:
                if 0 < self.stackslider.value() < self.imstack.nfiles:
                    self.stackslider.setValue(self.stackslider.value()-1)
        event.accept()


class TableView(QtWidgets.QTableWidget):
    def __init__(self, tablewidget, *args):
        QtWidgets.QTableWidget.__init__(self, *args)
        self.tablewidget = tablewidget
        self.data = []
        self.headerlabels = []
        self.setData()
        
    def setData(self):
        self.headerlabels = ['#']
        self.tablewidget.setColumnCount(1)
        self.tablewidget.setHorizontalHeaderLabels(self.headerlabels)
        for m, item in enumerate(self.data):
            row = self.tablewidget.rowCount()
            self.tablewidget.setRowCount(row+1)
            col = 0
            for el in item:
                cell = QtWidgets.QTableWidgetItem(str(el))
                self.tablewidget.setItem(row, col, cell)
                col += 1
        self.tablewidget.show()
        
    def addData(self, data, headerlabel):
        if headerlabel in self.headerlabels:
            pass  # deal with adding data to the same column
        self.headerlabels.append(headerlabel)
        col = self.tablewidget.columnCount()
        self.tablewidget.setColumnCount(col+1)
        self.tablewidget.setHorizontalHeaderLabels(self.headerlabels)
        if len(data) > self.tablewidget.rowCount():
            self.tablewidget.setRowCount(len(data))
        row = 0
        for el in data:
            cell = QtWidgets.QTableWidgetItem(str(round(el, 3)))
            self.tablewidget.setItem(row, col, cell)
            row += 1


class Pivdata:
    def __init__(self):
        self.x = []
        self.y = []
        self.u = []
        self.v = []
        self.ru = []
        self.rv = []
        self.width = 0
        self.height = 0
        self.nfields = 0
        self.centerx = 0
        self.centery = 0
        self.pixelsize = 1
        self.timeinterval = 1
        self.windowsize = 32
        self.halfwin = self.windowsize/2

    def rotacity(self, x, y, u, v, cx, cy):
        x = x-cx
        y = y-cy
        mag = np.linalg.norm(np.array([u, v]))
        theta = np.rad2deg(np.arctan2(y, x))+180
        fromN = np.deg2rad(270-theta)
        rotmat = np.array([[np.cos(fromN), -np.sin(fromN)], [np.sin(fromN), np.cos(fromN)]])
        rotuv = rotmat@np.array([[u], [v]])
        angle = np.rad2deg(np.arctan2(rotuv[1], rotuv[0])) + 180
        if angle < 90:
            xcomponent = -1 * (1 - angle / 90)
            ycomponent = -np.sqrt(1 - xcomponent ** 2)
        elif 90 <= angle < 180:
            xcomponent = (angle - 90) / 90
            ycomponent = -1 * np.sqrt(1 - xcomponent ** 2)
        elif 180 <= angle < 270:
            xcomponent = 1 - (angle - 180) / 90
            ycomponent = np.sqrt(1 - xcomponent ** 2)
        else:
            xcomponent = -1 * ((angle - 270) / 90)
            ycomponent = np.sqrt(1 - xcomponent ** 2)
        return xcomponent*mag, ycomponent*mag

    def lineariseField(self):
        self.ru = np.zeros((self.width, self.height, self.nfields))
        self.rv = np.zeros((self.width, self.height, self.nfields))
        for ff in range(self.nfields):
            for xx in range(self.height):
                for yy in range(self.width):
                    xcomp, ycomp = self.rotacity(self.x[xx, yy, ff]/self.halfwin, self.y[xx, yy, ff]/self.halfwin,
                                                 self.u[xx, yy, ff], -self.v[xx, yy, ff], self.centerx, self.centery)
                    self.ru[xx, yy, ff] = xcomp
                    self.rv[xx, yy, ff] = ycomp


def main():
    # Initialise GUI
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
        app.setWindowIcon(QIcon('PyV/icons/logo-02.png'))
    else:
        app = QtWidgets.QApplication.instance()
    app.setQuitOnLastWindowClosed(True)

    mainwindow = MainWindow()
    mainwindow.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
