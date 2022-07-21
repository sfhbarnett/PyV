from PyV.PIV import PIV
from tiffstack import TiffStack
import matplotlib.pyplot as plt
from piv_filters import localfilt
from naninterp import naninterp
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt6 import QtCore, QtWidgets
import sys
from superqt import QLabeledRangeSlider

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("PyV")
        self.layout = QtWidgets.QHBoxLayout()
        self.leftcontrols = QtWidgets.QVBoxLayout()
        self.centralcanvas = QtWidgets.QVBoxLayout()
        self.rightcontrols = QtWidgets.QVBoxLayout()

        # Left controls
        self.toggleroi = QtWidgets.QPushButton('Roi OFF ')
        #self.toggleroi.clicked.connect(self.toggleroimode)
        self.toggleroi.setToolTip('Toggle ROI mode on/off')
        self.roimode = 0
        self.driftcorbutton = QtWidgets.QPushButton('Drift')
        #self.driftcorbutton.clicked.connect(self.correctdrift)
        self.driftcorbutton.setToolTip('Correct the drift based on manually introduced point ROIs throughout the stack')
        self.PCCbutton = QtWidgets.QPushButton('PCC')
        #self.PCCbutton.clicked.connect(self.pccbuttonfunction)
        self.PCCbutton.setToolTip('Applies subpixel phase cross correlation to estimate drift')
        self.driftcheckbox = QtWidgets.QCheckBox("Apply drift", self)
        self.driftcheckbox.setToolTip("If a drift estimation has been made, this toggles the correction to the displayed"
                                      " data")
        self.driftcheckbox.setEnabled(False)
        #self.driftcheckbox.clicked.connect(self.viewdrift)
        self.autocontrast = QtWidgets.QCheckBox("AutoContrast", self)
        self.autocontrast.setToolTip("Toggle autocontrast on/off")
        self.autocontrast.setChecked(True)

        self.buttonbox = QtWidgets.QVBoxLayout()
        self.buttonbox.addStretch(1)
        self.buttonbox.addWidget(self.toggleroi)
        self.buttonbox.addWidget(self.driftcorbutton)
        self.buttonbox.addWidget(self.PCCbutton)
        self.buttonbox.addWidget(self.driftcheckbox)
        self.buttonbox.addWidget(self.autocontrast)
        self.buttonbox.addStretch(1)
        self.leftcontrols.addLayout(self.buttonbox)

        # Central Image controls

        self.sc = MplCanvas(self, width=7, height=8, dpi=100)
        self.filename = None
        path = '/Users/sbarnett/Documents/PIVData/fatima/ForSam/monolayer2/C1-20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi - 20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi #21.tif'
        self.imstack = TiffStack(path)
        self.plothandle = None
        self.sc.axes.get_xaxis().set_visible(False)
        self.sc.axes.get_yaxis().set_visible(False)
        self.sc.axes.set_aspect('auto')
        self.sc.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.s = self.sc.axes.scatter([], [], facecolors='none', edgecolors='w')
        self.mpl_toolbar = NavigationToolbar2QT(self.sc, None)
        self.contrastslider = QLabeledRangeSlider(QtCore.Qt.Vertical)
        self.contrastslider.setHandleLabelPosition(QLabeledRangeSlider.LabelPosition.LabelsBelow)
        self.contrastslider.label_shift_x = 10
        self.contrastslider.setRange(0, 200)
        #self.contrastslider.valueChanged.connect(self.update_contrast)
        self.imagecontrols = QtWidgets.QHBoxLayout()
        self.imagecontrols.addWidget(self.contrastslider)
        self.imagecontrols.addWidget(self.sc)

        self.centralcanvas.addLayout(self.imagecontrols)

        # Right controls

        # Combine
        self.layout.addLayout(self.leftcontrols)
        self.layout.addLayout(self.centralcanvas)
        self.layout.addLayout(self.rightcontrols)
        widget = QtWidgets.QWidget()
        widget.setLayout(self.layout)
        self.setCentralWidget(widget)
        self.setGeometry(100, 100, 1200, 900)




path = '/Users/sbarnett/Documents/PIVData/fatima/ForSam/monolayer2/C1-20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi - 20210708_MCF10ARAB5A_H2BGFP_Monolayer_Doxy_withoutDoxy.czi #21.tif'

tf = TiffStack(path)

# plt.imshow(tf.getimage(0))
# plt.show()
overlap = 0.5
windowsize = 32


# x, y, u, v = PIV(tf.getimage(0), tf.getimage(1), windowsize, overlap)
#
# u, v = localfilt(x, y, u, v, 2)
# u = naninterp(u)
# v = naninterp(v)
#
# plt.cla()
# plt.quiver(x, y, u, v, scale_units='dots', width=0.001, headlength=5, headwidth=3)
# plt.pause(0.1)
# plt.axis('equal')
# plt.show()

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