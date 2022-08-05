import os

from PyV.PIV import PIV
from tiffstack import tiffstack
import matplotlib.pyplot as plt
from piv_filters import localfilt
from naninterp import naninterp
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtGui import QIcon, QIntValidator, QDoubleValidator
import sys
import numpy as np
from main_gui import Ui_MainWindow
from superqt import QLabeledDoubleRangeSlider
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
        return self.x, self.y, self.u, self.v


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
        self.tableview = TableView(self.table)
        self.colormaps = ['jet','viridis','plasma','inferno']
        self.solids = ['black','green','red','cyan','magenta','yellow']
        self.solidcode = {'black':'k','green':'g','red':'r','cyan':'c','magenta':'m','yellow':'y'}
        self.fieldcolormapinput.addItems(self.colormaps+self.solids)
        self.quivercmap = 'jet'
        self.quivercolor = None
        self.progressbar = QtWidgets.QProgressBar()
        self.statusbar.addPermanentWidget(self.progressbar)
        self.alignmentbutton.setStyleSheet("")
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
        self.timeintervalinput.editingFinished.connect(self.settimeinterval)
        self.vrmsbutton.clicked.connect(self.calculatevrms)
        self.linearorderbutton.clicked.connect(self.calculateLinearOrderParameter)
        self.fieldcolormapinput.currentIndexChanged.connect(self.changecmap)
        self.showimagecheck.stateChanged.connect(self.showImage)

    def get_file(self):
        self.mplwidget.canvas.axes.cla()
        self.filename = None
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', directory='~/Documents')
        if self.filename[0] != "":
            self.filename = self.filename[0]
            self.imstack = tiffstack(self.filename)
            self.extent = [0, self.imstack.width,0,self.imstack.height]
            self.imagehandle = self.mplwidget.canvas.axes.imshow(self.imstack.getimage(0),extent=self.extent)
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
        self.piv.x = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.y = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.u = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.piv.v = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        # self.piv = externalPIV(self.x, self.y, self.u, self.v, self.windowsize, overlap, self.imstack)
        # self.piv.progressstatus.connect(self.updateprogress)
        # self.piv.start()
        start = time.time()
        for frame in range(self.imstack.nfiles-1):
            x, y, u, v = PIV(self.imstack.getimage(frame), self.imstack.getimage(frame+1), self.windowsize, overlap)
            u, v = localfilt(x, y, u, v, 2)
            u = naninterp(u)
            v = naninterp(v)
            self.piv.x[:, :, frame] = x * self.pixelsize
            self.piv.y[:, :, frame] = y * self.pixelsize
            self.piv.u[:, :, frame] = u * self.pixelsize / self.timeinterval
            self.piv.v[:, :, frame] = v * self.pixelsize / self.timeinterval
        self.showQuiver = True
        print(time.time()-start)
        self.makeQuiver()


    def makeQuiver(self):
        if self.quiver is not None:
            self.quiver.remove()
        if self.currentimage != 0 and self.showQuiver == True:
            frame = self.currentimage-1
            M = np.sqrt(self.piv.u[:, :, frame]*self.piv.u[:, :, frame]+self.piv.v[:, :, frame]*self.piv.v[:, :, frame])
            clipM = np.clip(M, 0, self.maxdispspeed)
            if self.quivercmap is not None:
                self.quiver = self.mplwidget.canvas.axes.quiver(self.piv.x[:, :, frame]+1,
                                                                self.piv.y[:, :, frame]+1,
                                                                self.piv.u[:, :, frame]*self.arrowscale,
                                                                self.piv.v[:, :, frame]*self.arrowscale, clipM,
                                                                scale_units='xy', scale=1,
                                                                cmap=self.quivercmap)
            else:
                self.quiver = self.mplwidget.canvas.axes.quiver(self.piv.x[:, :, frame]+1,
                                                                self.piv.y[:, :, frame]+1,
                                                                self.piv.u[:, :, frame]*self.arrowscale,
                                                                self.piv.v[:, :, frame]*self.arrowscale,
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

    def setwindowsize(self):
        self.windowsize = int(self.windowsizeinput.text())
        self.setFocus()

    def setarrowscale(self):
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

    def changecmap(self, index):
        cmap = self.fieldcolormapinput.currentText()
        if cmap in self.colormaps:
            self.quivercmap = cmap
            self.quivercolor = None
        else:
            self.quivercmap = None
            self.quivercolor = cmap
        self.makeQuiver()

    def exportfields(self):
        path = '/Users/sbarnett/PycharmProjects/PyV/exportloc'
        fmt = '%d', '%d', '%1.3f', '%1.3f'
        for frame in range(self.piv.x.shape[2]):
            array = np.vstack((self.piv.x[:, :, frame].flatten('F'), self.piv.y[:, :, frame].flatten('F'),
                               self.piv.u[:, :, frame].flatten('F'), self.piv.v[:, :, frame].flatten('F'))).T
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
        if self.currentimage != 0:
            frame = self.currentimage-1
            averageU = np.mean(self.piv.u[:, :, frame])
            averageV = np.mean(self.piv.v[:, :, frame])
            meanvector = np.array([averageU,averageV])
            upper = np.vstack((self.piv.u[:, :, frame].flatten('F'),self.piv.v[:, :, frame].flatten('F'))).T@meanvector
            lower = np.sqrt(self.piv.u[:, :, frame].flatten('F')**2+self.piv.v[:, :, frame].flatten('F')**2) * np.sqrt(meanvector[0]**2+meanvector[1]**2)
            value = np.divide(upper, lower)
            map = np.reshape(value, (63, 63))
            map = resize(map, (self.piv.u.shape[0]*self.windowsize//2, self.piv.u.shape[1]*self.windowsize//2), order=0, preserve_range=True)
            self.imagehandle.set_data(map)
            self.contrastslider.setRange(-1, 1)
            self.contrastslider.setValue((-1, 1))
            self.contrastslider.setSingleStep(0.01)
            self.contrastslider.setDecimals(2)
            self.contrastslider.label_shift_x = 10
            self.extent = [self.windowsize/4, self.imstack.width-self.windowsize/4, self.windowsize/4, self.imstack.width-self.windowsize/4]
            self.imagehandle.set_extent(self.extent)
            self.imagehandle.set_cmap('viridis')
            self.imagehandle.set_clim(-1, 1)
            self.mplwidget.canvas.fig.canvas.draw()

    def orientation(self):
        if self.currentimage != 0:
            frame = self.currentimage-1
            rho = np.sqrt(self.piv.u[:, :, frame]**2+self.piv.v[:, :, frame]**2)
            phi = np.arctan2(-self.piv.v[:, :, frame], self.piv.u[:, :, frame])
            theta = resize(np.rad2deg(phi), (self.piv.u.shape[0]*self.windowsize//2, self.piv.u.shape[1]*self.windowsize//2), order=0, preserve_range=True)
            self.imagehandle.set_data(theta+180)
            self.imagehandle.set_clim(0, 360)
            self.extent = [self.windowsize/4, self.imstack.width-self.windowsize/4, self.windowsize/4, self.imstack.width-self.windowsize/4]
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
                array = np.loadtxt(os.path.join(path,file), delimiter=',')
                self.piv.x[:, :, counter] = np.reshape(array[:,0],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.piv.y[:, :, counter] = np.reshape(array[:,1],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.piv.u[:, :, counter] = np.reshape(array[:,2],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                self.piv.v[:, :, counter] = np.reshape(array[:,3],[int(np.max(array[:,0]/(array[0,0]+1)))+1,int(np.max(array[:,1]/(array[0,1]+1)))+1])
                counter += 1
        self.showQuiver = True
        self.makeQuiver()

    def calculatevrms(self):
        """Calculate root mean square velocity of the vectorfield"""
        self.vrms = []
        for frame in range(self.piv.u.shape[2]):
            velocities = np.sqrt(self.piv.u[:, :, frame]**2+self.piv.v[:, :, frame]**2)
            self.vrms.append(np.sqrt(np.nanmean(velocities**2)))
        self.tableview.addData(self.vrms,'Vrms')

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
        self.mplwidget.canvas.figure.canvas.draw_idle()

    def calculateLinearOrderParameter(self):
        self.LOP = []
        for frame in range(self.piv.u.shape[2]):
            mean_u = np.nanmean(self.piv.u[:, :, frame])
            mean_v = np.nanmean(self.piv.v[:, :, frame])
            mean_uv2 = np.nanmean(self.piv.u[:, :, frame]**2 + self.piv.v[:, :, frame]**2)
            self.LOP.append((mean_u**2+mean_v**2)/mean_uv2)
        self.tableview.addData(self.LOP, 'LinearOrder')

    def natsort(self,tosort):
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
            pass # deal with adding data to the same column
        self.headerlabels.append(headerlabel)
        col = self.tablewidget.columnCount()
        self.tablewidget.setColumnCount(col+1)
        self.tablewidget.setHorizontalHeaderLabels(self.headerlabels)
        if len(data) > self.tablewidget.rowCount():
            self.tablewidget.setRowCount(len(data))
        row = 0
        for el in data:
            cell = QtWidgets.QTableWidgetItem(str(round(el, 4)))
            self.tablewidget.setItem(row, col, cell)
            row += 1


class Pivdata():
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
        self.pixelsize = 1
        self.timeinterval = 1
        self.windowsize = 32

    def rotaticity(self, x, y, u, v, cx, cy):
        x = x-cx
        y = y-cy
        theta = np.rad2deg(np.arctan2(y,x))+180
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
        elif angle >= 270:
            xcomponent = -1 * ((angle - 270) / 90)
            ycomponent = np.sqrt(1 - xcomponent ** 2)
        return xcomponent, ycomponent

    def lineariseField(self):
        for ff in range(self.nfields):
            for xx in range(self.height):
                for yy in range(self.width):
                    xcomp, ycomp = self.rotacity(self.x[xx,yy,ff],self.y[xx,yy,ff],self.u[xx,yy,ff],self.v[xx,yy,ff],self.width/2,self.height/2)
                    self.ru[xx,yy,ff] = xcomp
                    self.rv[xx,yy,ff] = ycomp





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