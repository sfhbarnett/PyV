from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PyQt6 import QtWidgets
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class MplWidget(QtWidgets.QWidget):
    def __init__(self,parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.canvas.axes.get_xaxis().set_visible(False)
        self.canvas.axes.get_yaxis().set_visible(False)
        self.canvas.axes.set_facecolor('#F0F0F0')
        self.canvas.fig.patch.set_facecolor('#EFEFEF')
        self.canvas.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.canvas.axes.set_aspect('auto')
        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)