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


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("PyV")

        self.leftcontrols = QtWidgets.QVBoxLayout()
        self.centralcanvas = QtWidgets.QVBoxLayout()
        self.rightcontrols = QtWidgets.QVBoxLayout()

        self.mincontrast = 0
        self.maxcontrast = 64000

        # Left controls =========================
        self.flo = QtWidgets.QFormLayout()
        self.pixelsizeinput = QtWidgets.QLineEdit()
        self.pixelsizeinput.editingFinished.connect(self.setPixelSize)
        self.pixelsize = 1
        self.maxspeedinput = QtWidgets.QLineEdit()
        self.maxspeedinput.editingFinished.connect(self.updatedisplayspeed)
        self.maxdispspeed = 100
        self.flo.addRow("Pixel Size", self.pixelsizeinput)
        self.flo.addRow("Max Speed", self.maxspeedinput)
        self.runPIVbutton = QtWidgets.QPushButton('Run. PIV analysis')
        self.runPIVbutton.clicked.connect(self.runPIV)
        self.runPIVbutton.setToolTip('Perform PIV analysis')
        self.autocontrast = QtWidgets.QCheckBox("AutoContrast", self)
        self.autocontrast.setToolTip("Toggle autocontrast on/off")
        self.autocontrast.setChecked(True)

        self.buttonbox = QtWidgets.QVBoxLayout()
        self.buttonbox.addStretch(1)
        self.buttonbox.addLayout(self.flo)
        self.buttonbox.addWidget(self.runPIVbutton)
        self.buttonbox.addWidget(self.autocontrast)
        self.buttonbox.addStretch(1)
        self.leftcontrols.addLayout(self.buttonbox)

        # Central Image controls ==========================

        self.sc = MplCanvas(self, width=7, height=8, dpi=100)
        self.filename = None
        path = '/Users/sbarnett/Documents/PIVData/fatima/ForSam/monolayer2/' \
               'C1-20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi -' \
               ' 20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi #21.tif'
        self.imstack = tiffstack(path)
        self.plothandle = None
        self.sc.axes.get_xaxis().set_visible(False)
        self.sc.axes.get_yaxis().set_visible(False)
        self.sc.axes.set_aspect('auto')
        self.sc.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.s = self.sc.axes.scatter([], [], facecolors='none', edgecolors='w')
        self.quiver = None
        self.mpl_toolbar = NavigationToolbar2QT(self.sc, None)
        self.contrastslider = QLabeledRangeSlider(QtCore.Qt.Vertical)
        self.contrastslider.setHandleLabelPosition(QLabeledRangeSlider.LabelPosition.LabelsBelow)
        self.contrastslider.label_shift_x = 10
        self.contrastslider.setRange(0, 200)
        self.contrastslider.valueChanged.connect(self.update_contrast)
        self.contrastslider.setMinimumWidth(100)
        self.contrastslider.label_shift_x = 10
        self.contrastslider.setEdgeLabelMode(self.contrastslider.EdgeLabelMode.NoLabel)
        self.imagecontrols = QtWidgets.QHBoxLayout()
        self.imagecontrols.addWidget(self.contrastslider)
        self.imagecontrols.addWidget(self.sc)

        self.centralcanvas.addLayout(self.imagecontrols)

        # Right controls ======================

        # Toolbar===================
        self.mpl_toolbar = NavigationToolbar2QT(self.sc, None)
        zoomact = QAction(QIcon("PyV/icons/zoom2.png"), "zoom", self)
        zoomact.setCheckable(True)
        zoomact.setShortcut("Ctrl+z")
        zoomact.triggered.connect(self.mpl_toolbar.zoom)

        panact = QAction(QIcon("PyV/icons/pan.png"), "pan", self)
        panact.setCheckable(True)
        panact.setShortcut("Ctrl+p")
        panact.setChecked(False)
        panact.triggered.connect(self.mpl_toolbar.pan)

        openact = QAction(QIcon("PyV/icons/open.png"), "open file", self)
        openact.setCheckable(False)
        openact.triggered.connect(self.get_file)

        saveact = QAction(QIcon("PyV/icons/save.png"), "Save drift corrected data", self)
        saveact.triggered.connect(self.savedrift)

        homeact = QAction(QIcon("PyV/icons/home.png"), "Reset axes to original view", self)
        homeact.setCheckable(False)
        homeact.triggered.connect(self.mpl_toolbar.home)

        self.toolbar = self.addToolBar("zoom")
        self.toolbar.addAction(panact)
        self.toolbar.addAction(zoomact)
        self.toolbar.addAction(openact)
        self.toolbar.addAction(saveact)
        self.toolbar.addAction(homeact)

        # Stack slider ================

        sliderholder = QtWidgets.QHBoxLayout()
        self.stackslider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.stackslider.setRange(0, self.imstack.nfiles - 1)
        self.stackslider.valueChanged.connect(self.move_through_stack)
        self.currentimage = 0
        sliderholder.addWidget(self.stackslider)
        self.label = QtWidgets.QLabel('0', self)
        self.label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.label.setMinimumWidth(40)

        sliderholder.addWidget(self.label)

        # Combine =====================

        self.layout = QtWidgets.QVBoxLayout()
        self.centrallayout = QtWidgets.QHBoxLayout()
        self.centrallayout.addLayout(self.leftcontrols)
        self.centrallayout.addLayout(self.centralcanvas)
        self.centrallayout.addLayout(self.rightcontrols)
        self.layout.addWidget(self.toolbar)
        self.layout.addLayout(self.centrallayout)
        self.layout.addLayout(sliderholder)
        widget = QtWidgets.QWidget()
        widget.setLayout(self.layout)
        self.setCentralWidget(widget)
        self.setGeometry(100, 100, 1200, 900)

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
        windowsize = 32
        overlap = 0.5
        self.x = np.zeros((int(self.imstack.width // (windowsize*overlap)-1), int(self.imstack.width // (windowsize*overlap)-1), self.imstack.nfiles-1))
        self.y = np.zeros((int(self.imstack.width // (windowsize*overlap)-1), int(self.imstack.width // (windowsize*overlap)-1), self.imstack.nfiles-1))
        self.u = np.zeros((int(self.imstack.width // (windowsize*overlap)-1), int(self.imstack.width // (windowsize*overlap)-1), self.imstack.nfiles-1))
        self.v = np.zeros((int(self.imstack.width // (windowsize*overlap)-1), int(self.imstack.width // (windowsize*overlap)-1), self.imstack.nfiles-1))
        for frame in range(self.imstack.nfiles-1):
            x, y, u, v = PIV(self.imstack.getimage(frame), self.imstack.getimage(frame+1), windowsize, overlap)
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
                                              self.u[:, :, frame], self.v[:, :, frame], clipM,scale_units='xy',scale=1, cmap=plt.cm.jet)
        else:
            self.quiver = self.sc.axes.quiver([], [], [], [])
        self.sc.fig.canvas.draw()

    def updatedisplayspeed(self):
        value = float(self.maxspeedinput.text())
        self.maxdispspeed = value
        self.makeQuiver()

    def setPixelSize(self):
        value = float(self.pixelsizeinput.text())
        self.u /= self.pixelsize
        self.v /= self.pixelsize
        self.pixelsize = value
        self.u *= self.pixelsize
        self.v *= self.pixelsize
        self.makeQuiver()





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