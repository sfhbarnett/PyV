# Form implementation generated from reading ui file 'main_gui.ui'
#
# Created by: PyQt6 UI code generator 6.1.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1405, 847)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.tabs = QtWidgets.QTabWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabs.sizePolicy().hasHeightForWidth())
        self.tabs.setSizePolicy(sizePolicy)
        self.tabs.setMinimumSize(QtCore.QSize(250, 300))
        self.tabs.setObjectName("tabs")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.formLayoutWidget = QtWidgets.QWidget(self.tab)
        self.formLayoutWidget.setGeometry(QtCore.QRect(30, 20, 164, 193))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.windowsizeinput = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.windowsizeinput.setObjectName("windowsizeinput")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.ItemRole.FieldRole, self.windowsizeinput)
        self.label_3 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label_3)
        self.maxspeedinput = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.maxspeedinput.setObjectName("maxspeedinput")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.ItemRole.FieldRole, self.maxspeedinput)
        self.label_2 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label_2)
        self.pixelsizeinput = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.pixelsizeinput.setObjectName("pixelsizeinput")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.ItemRole.FieldRole, self.pixelsizeinput)
        self.label = QtWidgets.QLabel(self.formLayoutWidget)
        self.label.setObjectName("label")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label)
        self.label_4 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_4.setObjectName("label_4")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label_4)
        self.autocontrast = QtWidgets.QCheckBox(self.formLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.autocontrast.sizePolicy().hasHeightForWidth())
        self.autocontrast.setSizePolicy(sizePolicy)
        self.autocontrast.setText("")
        self.autocontrast.setChecked(True)
        self.autocontrast.setObjectName("autocontrast")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.ItemRole.FieldRole, self.autocontrast)
        self.arrowscaleinput = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.arrowscaleinput.setObjectName("arrowscaleinput")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.ItemRole.FieldRole, self.arrowscaleinput)
        self.label_5 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label_5)
        self.label_6 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_6.setObjectName("label_6")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.ItemRole.LabelRole, self.label_6)
        self.timeintervalinput = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.timeintervalinput.setObjectName("timeintervalinput")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.ItemRole.FieldRole, self.timeintervalinput)
        self.runPIVbutton = QtWidgets.QPushButton(self.tab)
        self.runPIVbutton.setGeometry(QtCore.QRect(50, 220, 113, 32))
        self.runPIVbutton.setObjectName("runPIVbutton")
        self.tabs.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.vrmsbutton = QtWidgets.QPushButton(self.tab_2)
        self.vrmsbutton.setGeometry(QtCore.QRect(10, 10, 104, 32))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.vrmsbutton.sizePolicy().hasHeightForWidth())
        self.vrmsbutton.setSizePolicy(sizePolicy)
        self.vrmsbutton.setMinimumSize(QtCore.QSize(104, 0))
        self.vrmsbutton.setObjectName("vrmsbutton")
        self.linearorderbutton = QtWidgets.QPushButton(self.tab_2)
        self.linearorderbutton.setGeometry(QtCore.QRect(120, 10, 104, 32))
        self.linearorderbutton.setMaximumSize(QtCore.QSize(104, 16777215))
        self.linearorderbutton.setObjectName("linearorderbutton")
        self.orientationbutton = QtWidgets.QPushButton(self.tab_2)
        self.orientationbutton.setGeometry(QtCore.QRect(120, 40, 104, 32))
        self.orientationbutton.setMinimumSize(QtCore.QSize(104, 0))
        self.orientationbutton.setMaximumSize(QtCore.QSize(104, 16777215))
        self.orientationbutton.setObjectName("orientationbutton")
        self.alignmentbutton = QtWidgets.QPushButton(self.tab_2)
        self.alignmentbutton.setGeometry(QtCore.QRect(10, 40, 104, 32))
        self.alignmentbutton.setObjectName("alignmentbutton")
        self.fieldcolormapinput = QtWidgets.QComboBox(self.tab_2)
        self.fieldcolormapinput.setGeometry(QtCore.QRect(120, 80, 104, 26))
        self.fieldcolormapinput.setObjectName("fieldcolormapinput")
        self.label_7 = QtWidgets.QLabel(self.tab_2)
        self.label_7.setGeometry(QtCore.QRect(15, 75, 92, 31))
        self.label_7.setObjectName("label_7")
        self.tabs.addTab(self.tab_2, "")
        self.horizontalLayout_2.addWidget(self.tabs)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.contrastslider = QLabeledRangeSlider(self.centralwidget)
        self.contrastslider.setOrientation(QtCore.Qt.Orientation.Vertical)
        self.contrastslider.setObjectName("contrastslider")
        self.horizontalLayout.addWidget(self.contrastslider)
        spacerItem1 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.mplwidget = MplWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplwidget.sizePolicy().hasHeightForWidth())
        self.mplwidget.setSizePolicy(sizePolicy)
        self.mplwidget.setMinimumSize(QtCore.QSize(400, 400))
        self.mplwidget.setObjectName("mplwidget")
        self.horizontalLayout.addWidget(self.mplwidget)
        self.table = QtWidgets.QTableWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.table.sizePolicy().hasHeightForWidth())
        self.table.setSizePolicy(sizePolicy)
        self.table.setMaximumSize(QtCore.QSize(300, 16777215))
        self.table.setObjectName("table")
        self.table.setColumnCount(0)
        self.table.setRowCount(0)
        self.horizontalLayout.addWidget(self.table)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.stackslider = QtWidgets.QSlider(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.stackslider.sizePolicy().hasHeightForWidth())
        self.stackslider.setSizePolicy(sizePolicy)
        self.stackslider.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.stackslider.setObjectName("stackslider")
        self.horizontalLayout_3.addWidget(self.stackslider)
        self.stackpos = QtWidgets.QLabel(self.centralwidget)
        self.stackpos.setObjectName("stackpos")
        self.horizontalLayout_3.addWidget(self.stackpos)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1405, 24))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(MainWindow)
        self.toolBar.setObjectName("toolBar")
        MainWindow.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.toolBar)
        self.actionzoom = QtGui.QAction(MainWindow)
        self.actionzoom.setCheckable(True)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/icons/zoom2.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionzoom.setIcon(icon)
        self.actionzoom.setObjectName("actionzoom")
        self.actionPan = QtGui.QAction(MainWindow)
        self.actionPan.setCheckable(True)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/icons/pan.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionPan.setIcon(icon1)
        self.actionPan.setObjectName("actionPan")
        self.actionopen = QtGui.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/icons/open.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionopen.setIcon(icon2)
        self.actionopen.setObjectName("actionopen")
        self.actionhome = QtGui.QAction(MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/icons/home.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionhome.setIcon(icon3)
        self.actionhome.setObjectName("actionhome")
        self.actionSave = QtGui.QAction(MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/icons/save.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionSave.setIcon(icon4)
        self.actionSave.setObjectName("actionSave")
        self.actionexport = QtGui.QAction(MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/icons/export.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionexport.setIcon(icon5)
        self.actionexport.setObjectName("actionexport")
        self.actionimport = QtGui.QAction(MainWindow)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap(":/icons/import.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionimport.setIcon(icon6)
        self.actionimport.setObjectName("actionimport")
        self.menubar.addAction(self.menuFile.menuAction())
        self.toolBar.addAction(self.actionopen)
        self.toolBar.addAction(self.actionSave)
        self.toolBar.addAction(self.actionexport)
        self.toolBar.addAction(self.actionimport)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionzoom)
        self.toolBar.addAction(self.actionPan)
        self.toolBar.addAction(self.actionhome)

        self.retranslateUi(MainWindow)
        self.tabs.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.windowsizeinput.setToolTip(_translate("MainWindow", "PIV window size (pixels)"))
        self.windowsizeinput.setText(_translate("MainWindow", "32"))
        self.label_3.setText(_translate("MainWindow", "Max Speed"))
        self.maxspeedinput.setToolTip(_translate("MainWindow", "Quiver color will be clipped to this value"))
        self.maxspeedinput.setText(_translate("MainWindow", "100"))
        self.label_2.setText(_translate("MainWindow", "Pixel Size"))
        self.pixelsizeinput.setToolTip(_translate("MainWindow", "Pixel size of the data, PIV results will use this as the scale"))
        self.pixelsizeinput.setText(_translate("MainWindow", "1"))
        self.label.setText(_translate("MainWindow", "Window Size"))
        self.label_4.setText(_translate("MainWindow", "autocontrast"))
        self.arrowscaleinput.setToolTip(_translate("MainWindow", "Set arrow scaling for visualisation"))
        self.arrowscaleinput.setText(_translate("MainWindow", "1"))
        self.label_5.setText(_translate("MainWindow", "Arrow Scale"))
        self.label_6.setText(_translate("MainWindow", "Time interval"))
        self.timeintervalinput.setText(_translate("MainWindow", "1"))
        self.runPIVbutton.setToolTip(_translate("MainWindow", "Perform PIV analysis with current settings"))
        self.runPIVbutton.setText(_translate("MainWindow", "Run PIV"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab), _translate("MainWindow", "PIV"))
        self.vrmsbutton.setText(_translate("MainWindow", "Vrms"))
        self.linearorderbutton.setText(_translate("MainWindow", "Linear Order"))
        self.orientationbutton.setText(_translate("MainWindow", "Orientation"))
        self.alignmentbutton.setText(_translate("MainWindow", "Alignment"))
        self.label_7.setText(_translate("MainWindow", "Field Colormap"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_2), _translate("MainWindow", "Analysis"))
        self.stackpos.setText(_translate("MainWindow", "0"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar"))
        self.actionzoom.setText(_translate("MainWindow", "zoom"))
        self.actionzoom.setToolTip(_translate("MainWindow", "zooms on image"))
        self.actionPan.setText(_translate("MainWindow", "Pan"))
        self.actionopen.setText(_translate("MainWindow", "open"))
        self.actionhome.setText(_translate("MainWindow", "home"))
        self.actionSave.setText(_translate("MainWindow", "Save"))
        self.actionexport.setText(_translate("MainWindow", "export"))
        self.actionexport.setToolTip(_translate("MainWindow", "export quiver fields"))
        self.actionimport.setText(_translate("MainWindow", "import"))
        self.actionimport.setToolTip(_translate("MainWindow", "Import fields"))
from mplwidget import MplWidget
from superqt import QLabeledRangeSlider


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
