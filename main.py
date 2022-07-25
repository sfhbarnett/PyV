from PyV.PIV import PIV
from tiffstack import tiffstack
import matplotlib.pyplot as plt
from piv_filters import localfilt
from naninterp import naninterp
from matplotlib.figure import Figure
import matplotlib.cm as cm
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtGui import QAction, QIcon
import sys
from superqt import QLabeledRangeSlider
import numpy as np
from main_gui import Ui_MainWindow


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.seticons()
        self.actionopen.triggered.connect(self.get_file)

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

    def get_file(self):
        self.sc.axes.cla()
        self.filename = None
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', directory='~/Documents')
        if self.filename[0] != "":
            self.filename = self.filename[0]
            self.imstack = tiffstack(self.filename)
            self.plothandle = self.sc.axes.imshow(self.imstack.getimage(0))
            self.stackslider.setRange(0, self.imstack.nfiles - 1)
            self.stackslider.setValue(0)
            self.plothandle.set_cmap('gray')
            self.sc.fig.canvas.draw()
            self.contrastslider.setRange(0, self.imstack.maximum * 1.5)
            self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))

    def savedrift(self):
        pass

    def move_through_stack(self, value):
        """ Updates the current image in the viewport"""
        self.currentimage = value

        self.plothandle.set_data(self.imstack.getimage(value))

        # Check and update contrast
        if not self.autocontrast.isChecked():
            self.plothandle.set_clim([self.mincontrast, self.maxcontrast])
        else:
            self.contrastslider.setRange(0, self.imstack.maximum*1.5)
            self.contrastslider.setValue((self.imstack.minimum, self.imstack.maximum))
        if self.quiver is not None:
            self.makeQuiver()
        self.label.setText(str(value))

    def update_contrast(self, value):
        self.mincontrast = value[0]
        self.maxcontrast = value[1]
        if self.plothandle:
            self.plothandle.set_clim(self.mincontrast, self.maxcontrast)
            self.sc.fig.canvas.draw()

    def runPIV(self):
        overlap = 0.5
        self.x = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.y = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.u = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        self.v = np.zeros((int(self.imstack.width // (self.windowsize*overlap)-1), int(self.imstack.width // (self.windowsize*overlap)-1), self.imstack.nfiles-1))
        for frame in range(self.imstack.nfiles-1):
            x, y, u, v = PIV(self.imstack.getimage(frame), self.imstack.getimage(frame+1), self.windowsize, overlap)
            u, v = localfilt(x, y, u, v, 2)
            u = naninterp(u)
            v = naninterp(v)
            self.x[:, :, frame] = x
            self.y[:, :, frame] = y
            self.u[:, :, frame] = u
            self.v[:, :, frame] = v
        self.makeQuiver()

    def makeQuiver(self):
        if self.quiver is not None:
            self.quiver.remove()
        if self.currentimage != 0:
            frame = self.currentimage-1
            M = np.sqrt(self.u[:, :, frame]*self.u[:, :, frame]+self.u[:, :, frame]*self.v[:, :, frame])
            clipM = np.clip(M, 0, self.maxdispspeed)
            self.quiver = self.sc.axes.quiver(self.x[:, :, frame], self.y[:, :, frame],
                                              self.u[:, :, frame], self.v[:, :, frame], clipM,
                                              scale_units='xy',scale=1, cmap=plt.cm.jet)
        else:
            self.quiver = self.sc.axes.quiver([], [], [], [])
        self.sc.fig.canvas.draw()

    def updatedisplayspeed(self):
        value = float(self.maxspeedinput.text())
        self.maxdispspeed = value
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