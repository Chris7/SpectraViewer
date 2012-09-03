# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SpecView.ui'
#
# Created: Sun Sep 02 19:37:23 2012
#      by: PyQt4 UI code generator 4.9.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(737, 397)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.verticalLayout.addWidget(self.tabWidget)
        self.gridLayout.addLayout(self.verticalLayout, 1, 0, 1, 1)
        self.instructions = QtGui.QLabel(self.centralwidget)
        self.instructions.setMinimumSize(QtCore.QSize(0, 100))
        self.instructions.setTextFormat(QtCore.Qt.AutoText)
        self.instructions.setScaledContents(False)
        self.instructions.setAlignment(QtCore.Qt.AlignCenter)
        self.instructions.setObjectName(_fromUtf8("instructions"))
        self.gridLayout.addWidget(self.instructions, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 737, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(-1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Spectra Viewer", None, QtGui.QApplication.UnicodeUTF8))
        self.instructions.setText(QtGui.QApplication.translate("MainWindow", "To load a file, simply drag and drop it!", None, QtGui.QApplication.UnicodeUTF8))

