#!/usr/bin/env python

# embedding_in_qt4.py --- Simple Qt4 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#               2006 Darren Dale
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.

import sys, os, random
sys.path.append('/sw/lib/qt4-mac/lib/python2.7/site-packages/')

from PyQt4 import QtGui, QtCore

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from mutinf_analysis import *

progname = os.path.basename(sys.argv[0])
progversion = "0.1"


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None):
        fig = Figure()
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self):
        j = mutInfmat('2esk_demo.txt',[])
        j.unsort(self.axes) 
        

class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        #timer = QtCore.QTimer(self)
        #QtCore.QObject.connect(timer, QtCore.SIGNAL("timeout()"), self.update_figure)
        #timer.start(1000)

    def update_figure(self,pos1,pos2):

        self.axes.scatter(int(pos1),int(pos2))
        self.draw()


class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QGridLayout(self.main_widget)
        sc = MyStaticMplCanvas(self.main_widget) #, width=8, height=8)
        self.dc = MyDynamicMplCanvas(self.main_widget) #, width=3, height=10)
        
	self.edit1 = QtGui.QLineEdit(self)
	label1= QtGui.QLabel(self)
	label1.setText('Position 1')
	label1.setAlignment(QtCore.Qt.AlignRight)

	self.edit2 = QtGui.QLineEdit(self)
	label2= QtGui.QLabel(self)
	label2.setAlignment(QtCore.Qt.AlignRight)
	label2.setText('Position 2')
	
	l.addWidget(sc,1,1,2,1)
        l.addWidget(self.dc,1,1,1,2)

	l.addWidget(self.edit1,2,)
        l.addWidget(label1,1,0)
	l.addWidget(self.edit2,2,1)
	l.addWidget(label2,2,0)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("Just testing here.")
    
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Return:
	    print 'Position 1 is:', self.edit1.text()
            print 'Position 2 is:', self.edit2.text()
	    #if self.edit.text()

        self.dc.update_figure(self.edit1.text(),self.edit2.text())

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
