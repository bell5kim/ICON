import os, sys

# from PyQt4 import QtCore, QtGui
# Import PyQt Widgets and Matplotlib canvas for actually used PyQt version
if "PyQt5" in sys.modules:
    from PyQt5.QtWidgets import QSizePolicy, QWidget, QVBoxLayout, QAction, QApplication
    from PyQt5.QtCore import QSize, QObject, pyqtSignal
    from PyQt5.QtGui import QImage
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.backend_bases import NavigationToolbar2
elif "PyQt4" in sys.modules:
    from PyQt4.QtGui import QSizePolicy, QWidget, QVBoxLayout, QAction, QApplication, QImage
    from PyQt4.QtCore import QSize, QObject, pyqtSignal
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.backend_bases import NavigationToolbar2
else:
    raise SystemError("PyQt4 or PyQt5 need to be imported first")

from mpl_toolkits.mplot3d import Axes3D

class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, width = 10, height = 10, dpi = 100, sharex = None, sharey = None):
        self.fig = Figure(figsize = (width, height), dpi=dpi, facecolor = '#FFFFFF')
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95)
        self.xtitle="x"
        self.ytitle="y"
        self.PlotTitle = "TiFF Image"
        self.grid_status = True
        self.xaxis_style = 'linear'
        self.yaxis_style = 'linear'
        self.format_labels()
        self.ax.hold(True) 
            
        FigureCanvas.__init__(self, self.fig)
        #self.fc = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def format_labels(self):
        self.ax.set_title(self.PlotTitle)
        self.ax.title.set_fontsize(5)
        self.ax.set_xlabel(self.xtitle, fontsize = 5)
        self.ax.set_ylabel(self.ytitle, fontsize = 5)
        labels_x = self.ax.get_xticklabels()
        labels_y = self.ax.get_yticklabels()

        for xlabel in labels_x:
            xlabel.set_fontsize(5)
        for ylabel in labels_y:
            ylabel.set_fontsize(5)
            ylabel.set_color('b')

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(10, 10)


class QtMplCanvas(FigureCanvas):
    def __init__(self, parent=None, width = 6.5, height = 5.5, dpi = 100, sharex = None, sharey = None, fig = None):
        if fig == None:
            self.fig = Figure(figsize = (width, height), dpi=dpi, facecolor = '#FFFFFF')
            self.ax = self.fig.add_subplot(111, projection='3d')
            self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
            self.ax.hold(True)
        else:
            self.fig = fig

        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
            QSizePolicy.Expanding,
            QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def format_labels(self):
        self.ax.set_title(self.PlotTitle)
        self.ax.title.set_fontsize(5)
        self.ax.set_xlabel(self.xtitle, fontsize = 5)
        self.ax.set_ylabel(self.ytitle, fontsize = 5)
        self.ax.set_zlabel(self.ztitle, fontsize = 5)
        labels_x = self.ax.get_xticklabels()
        labels_y = self.ax.get_yticklabels()
        labels_z = self.ax.get_zticklabels()

        for xlabel in labels_x:
            xlabel.set_fontsize(5)
        for ylabel in labels_y:
            ylabel.set_fontsize(5)
            ylabel.set_color('b')  
        for zlabel in labels_z:
            zlabel.set_fontsize(5)        

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(10, 10)
    
        
class MPL_Widget(QWidget):
    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self.canvas = MyMplCanvas()
        # self.toolbar = NavigationToolbar(self.canvas, self.canvas)
        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.canvas)
        # self.vbox.addWidget(self.toolbar)
        self.setLayout(self.vbox)
        

class MPL_Widget_3D(QWidget):
    def __init__(self, parent = None, enableAutoScale = False, enableCSV = False, enableEdit = False, fig = None):
        QWidget.__init__(self, parent)
        self.canvas = QtMplCanvas(fig)
        self.canvas.ax.mouse_init()
        self.toolbar = NavigationToolbar(self.canvas, self.canvas)
        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.canvas)
        self.vbox.addWidget(self.toolbar)
        self.setLayout(self.vbox)

        ###########SAVING FIGURE TO CLIPBOARD##########
        self.cb = None #will be used for the clipboard
        self.tempPath = getHomeDir()
        self.tempPath = os.path.join(self.tempPath,'tempMPL.png')

        self.mpl2ClipAction = QAction("Save to Clipboard",  self)
        self.mpl2ClipAction.setShortcut("Ctrl+C")
        self.addAction(self.mpl2ClipAction)
        # QObject.connect(self.mpl2ClipAction,pyqtSignal("triggered()"), self.mpl2Clip)
        self.mpl2ClipAction.triggered.connect(self.mpl2Clip)

    def mpl2Clip(self):
        try:
            self.canvas.fig.savefig(self.tempPath)
            tempImg = QImage(self.tempPath)
            self.cb = QApplication.clipboard()
            self.cb.setImage(tempImg)
        except:
            print 'Error copying figure to clipboard'
            errorMsg = "Sorry: %s\n\n:%s\n"%(sys.exc_type, sys.exc_value)
            print errorMsg
            
            
def valid(path):
    if path and os.path.isdir(path):
        return True
    return False

def env(name):
    return os.environ.get( name, '' )

def getHomeDir():
    if sys.platform != 'win32':
        return os.path.expanduser( '~' )

    homeDir = env( 'USERPROFILE' )
    if not valid(homeDir):
        homeDir = env( 'HOME' )
        if not valid(homeDir) :
            homeDir = '%s%s' % (env('HOMEDRIVE'),env('HOMEPATH'))
            if not valid(homeDir) :
                homeDir = env( 'SYSTEMDRIVE' )
                if homeDir and (not homeDir.endswith('\\')) :
                    homeDir += '\\'
                if not valid(homeDir) :
                    homeDir = 'C:\\'
    return homeDir
