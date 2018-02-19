#!/usr/bin/env python

try:
    from PyQt4.uic import loadUiType
except:
    from PyQt5.uic import loadUiType
    
import sys, os, re

try:
    # from PyQt4 import QtCore, QtGui, QtXml
    from PyQt4.QtGui import QMainWindow, QSizePolicy, QWidget, QVBoxLayout, QAction, QApplication, QImage, QTableWidgetItem,\
                            QColor, QMessageBox, QFileDialog
    from PyQt4.QtCore import QSize, QObject, pyqtSignal, QTime, QStringList, QFile, QTextStream, QString
    from PyQt4.QtXml import QDomDocument
    print("Running Icon Gating with PyQt4...")
except:
    from PyQt5.QtWidgets import QMainWindow, QSizePolicy, QWidget, QVBoxLayout, QAction, QApplication, QTableWidgetItem, \
                            QMessageBox, QFileDialog
    from PyQt5.QtCore import QSize, QObject, pyqtSignal, QTime, QStringListModel, QFile, QTextStream
    from PyQt5.QtGui import QImage, QColor
    from PyQt5.QtXml import QDomDocument
    print("Running Icon Gating with PyQt5...")    
    
from matplotlib.patches import Ellipse

import numpy as NUMPY
import datetime as DATETIME


# from ui_BrachyPlan import Ui_MainWindow
try: 
    from ui_BrachyPlan import Ui_MainWindow
except:
    Ui_MainWindow, QMainWindow = loadUiType('IconGating.ui')

        
class Plot_Widget(QMainWindow, Ui_MainWindow):
    
    def __init__(self, data2plot=None, parent=None):
        super(Plot_Widget, self).__init__(parent)
        self.setupUi(self)
        
        self.initSettings()
        
        self.statusBar.showMessage('Ready')
        
        self.ExitPushButton.clicked.connect(self.exitWindow)
        self.HomePushButton.clicked.connect(self.setHome)  
        self.SavePushButton.clicked.connect(self.saveSettingsXml) 
        self.GKDirListWidget.itemClicked.connect(self.updateGKLogList)
        self.GKLogListWidget.itemClicked.connect(self.readGKLog)
        self.TableWidgetEvent.itemSelectionChanged.connect(self.setEventRange)
        
        self.comboBoxSort.activated.connect(self.sortEvent)
        self.checkBoxAspect.stateChanged.connect(self.updatePlot2dError)
        self.checkBoxStdRef.stateChanged.connect(self.updatePlot2dError)
        self.checkBoxStdMean.stateChanged.connect(self.updatePlot2dError)
        self.checkBoxAnnotate.stateChanged.connect(self.updatePlot2dError)
        self.checkBoxSimple.stateChanged.connect(self.simpleEventView)
        
        self.pushButtonXY.clicked.connect(self.rotate3dXY)
        self.pushButtonYZ.clicked.connect(self.rotate3dYZ)
        self.pushButtonXZ.clicked.connect(self.rotate3dXZ)
        
        self.ReportPushButton.clicked.connect(self.reportShots)     
        
        
        self.readSettingsXml()
        


            
    def exitWindow(self):
        self.close()


    def initSettings(self):

        self.ROOT = ''
        self.OS = sys.platform
        slash = '/'
        
        # XML Path Setting 
        if self.OS == 'linux2':
            self.HOME = os.environ["HOME"]
        elif self.OS == 'win32':
            self.HOME = os.environ["USERPROFILE"]
            slash = '\\'
        else:
            self.HOME = os.environ["HOMEDRIVE"] + os.environ["HOMEPATH"]
        
        print "self.HOME = ", self.HOME   
        self.initFile = self.HOME +slash+ "IconGating.xml"
        print "self.initFile = ", self.initFile  
        
        
        print 'Today = ', DATETIME.date.today()
        # print 'Today = ', DATETIME.datetime.today()
        # XML Setting
        self.domDocument = QDomDocument("IconGating")
        
        self.xMouse, self.yMouse, self.zMouse = 0, 0, 0
        self.xPos, self.yPos, self.zPos = 0, 0, 0
        
        self.patientDir = ""
        self.patient = ""

        self.xList = []
        self.yList = []
        self.dim = []
        self.RS = []
        
                
    def initTableShot(self):
                  
        # Remove All Rows in the Table
        nRows = self.TableWidgetShot.rowCount()
        for i in range(0, nRows):
            self.TableWidgetShot.removeRow(nRows - i - 1)
        
        # Remove All Columns in the TAble
        nColumns = self.TableWidgetShot.columnCount()
        for i in range(0, nColumns):
            self.TableWidgetShot.removeColumn(nColumns - i - 1)
        
        #if "PyQt5" in sys.modules:
        #    self.TableLabelShot = QStringListModel()
        #if "PyQt4" in sys.modules:
        #    self.TableLabelShot = QStringList() 

        #self.TableLabelShot << "ShotID" << "Planned" << "Elapsed" \
        #<< "Started" << "Stopped" << "Beam On" << "Beam Off" << "Duration" \
        #<< "R_Mean" << "R_StdDev" << "X_Mean" << "X_StdDev" \
        #<< "Y_Mean" << "Y_StdDev" << "Z_Mean" << "Z_StdDev"
        
        self.TableLabelShot = [
            "ShotID", "Planned", "Elapsed", "Started", "Stopped", "Beam On",
            "Beam Off", "Duration", "R_Mean", "R_StdDev", "X_Mean", "X_StdDev",
            "Y_Mean", "Y_StdDev", "Z_Mean", "Z_StdDev"]

        self.TableWidgetShot.setColumnCount(len(self.TableLabelShot))
        self.TableWidgetShot.setHorizontalHeaderLabels(self.TableLabelShot)    
  
            
                  
    def sortEvent(self):
        i = self.comboBoxSort.currentIndex()
        self.TableWidgetEvent.sortItems(i) 
        # print i, 'sorting'
        
        
    def MeanStdOverall (self, X, Y, M = 0.0):
                
        if M == 0.0:
            # Calculate Mean by Trapezoidal rule
            M = NUMPY.trapz(Y, x=X)/(max(X)-min(X))
            if X[0] > X[1]: print X
            if M < 0:
                # print X, Y
                print 'ERROR: Negative Value of M = ', M, NUMPY.trapz(Y, x=X)
        
        if len(Y) == 1:
            return[M, 0.0]
                        
        S2 = 0.0
        for i in range(len(X)-1):
            a, b = X[i], X[i+1]
            A = float(Y[i+1]-Y[i])/float(X[i+1]-X[i])
            B = Y[i]-A*X[i]-M
            
            # m = A/2*(b+a)+B
            if A == 0: 
                s2 = (b-a)*B*B
            else:
                s2 = (NUMPY.power(A*b+B,3.0) - NUMPY.power(A*a+B, 3.0))/(3*A)
            
            # print X[i], X[i+1], Y[i], Y[i+1]
            # print 'i =\t',i,'\tA =\t', A,'\tB =\t', B, '\tm =\t', m,'\ts2 =\t', s2

            S2 += s2
            if S2 < 0:
                print 'i =\t',i,'\tA =\t', A,'\tB =\t', B, '\tM =\t', M,'\ts2 =\t', s2
        
        # print 'M, S2, min(X), max(X) = ', M, S2, min(X), max(X)
        
        return [M, NUMPY.sqrt(S2/(max(X)-min(X)))]
    
    
    def MeanStdOverallx (self, T, X, M = 0.0):
        if M == 0.0:
            # Calculate Mean by Trapezoidal rule
            # print ('len(T) = ',len(T),'max(T) = ', max(T),' min(T) = ',min(T))
            M = NUMPY.trapz(X, x=T)/(max(T)-min(T))
            
        D = X - M
        # S = NUMPY.sqrt(NUMPY.trapz(D, x=T)/(max(T)-min(T)))
        
        # Absolute Value and Integrate
        Tabs, Xabs = [], []
        for i in range(len(D)-1):
            # print X[i], Ydev[i]
            Tabs.append(T[i])
            Xabs.append(D[i])
            if D[i]*D[i+1] < 0:
                m = (D[i+1]-D[i])/(T[i+1]-T[i])
                t = T[i] - D[i]/m
                # print x, Ydev[i]+m*(x-X[i])   
                Tabs.append(t)
                # Xabs.append(D[i]+m*(t-T[i]))
                Xabs.append(0.0)
        
        Tabs.append(T[len(T)-1])
        Xabs.append(D[len(X)-1])
            
        # print '*** MeanStdOverall:: len(T), Tabs = ', len(T), Tabs
        #f or i in range(len(Xabs)):
        #     print Xabs[i], Yabs[i]
            
        if len(T) == 1:
            return[M,0.0]
            

        # print '*** MeanStdOverall:: Min, Max = ', min(Tabs), max(Tabs), M+min(Xabs), M+max(Xabs)
        S = NUMPY.trapz(NUMPY.abs(Xabs), x=Tabs)/(max(Tabs)-min(Tabs))


        return [M, S]
            
            
    def time2sec (self, h, m, s, ms):
        return h*60*60 + m*60 + s + ms/1000.0
        
            
    def lineIntp(self, x, P1, P2):
        if P1[0] is P2[0]:
            return 10000.0
        
        return (P2[1]-P1[1])/(P2[0]-P1[0])*(x - P1[0]) + P1[1]          
    
    def MaxMinEachBeam(self, Xin, Yin):

        X = NUMPY.asarray(Xin)-self.T0
        Y = Yin    
        
        return [[min(X), max(X)], [min(Y), max(Y)]]
    
    def MeanStdEachBeam(self, Xin, Yin, M=0.0):
        # print " "
        # print "*** self.BeamStateValueList = ", self.BeamStateValueList
        # print "*** self.BeamStateTimeList  = ", self.BeamStateTimeList
        # print "*** self.T0 = ", self.T0
        # print "*** Xin = ", Xin
        # print "*** Yin = ", Yin
        X = NUMPY.asarray(Xin)-self.T0
        Y = Yin
        
        # print 'X = ', list(X)
        # print 'Y = ', Y
        # print "*** MeanStdEachBeam:: X[0] and max(X) = ", X[0], max(X)    
        
        BeamMeanStd = []
        
        nBeamState = len(self.BeamStateValueList)
        if nBeamState == 0: return[[],[],[]]
        # print '*** MeanStdEachBeam:: nBeamState = ', nBeamState
        BeamStateTimeList = NUMPY.asarray(self.BeamStateTimeList)-self.T0
        # print '*** MeanStdEachBeam:: len(BeamStateTimeList) = ', len(BeamStateTimeList)
        # print '*** MeanStdEachBeam:: BeamStateTimeList) = ', BeamStateTimeList

        BeamOnOffList = []
        BeamOn, BeamOff = 0.0, 0.0
        for i in range(nBeamState):
            # print 'BeamStateTime[',i,'] = ', BeamStateTimeList[i], self.BeamStateValueList[i]
            if self.BeamStateValueList[i].find('ON')  == 0:
                BeamOn = BeamStateTimeList[i]
                # print i, 'BeamOn = ', BeamOn
                
            if self.BeamStateValueList[i].find('OFF') == 0:
                BeamOff = BeamStateTimeList[i]
                # print i, 'BeamOff = ', BeamOff
                                        
                if BeamOff > BeamOn and BeamOn > X[0] and BeamOn < max(X): 
                    BeamOnOffList.append([BeamOn, BeamOff])
                    
        # print '***** MeanStdEachBeam:: len(BeamOnOffList) = ', len(BeamOnOffList)  
        # print '***** MeanStdEachBeam:: BeamOnOffList = ', BeamOnOffList          

        nX = len(X)                   # Number of X 
        Xa = list(X)
        Ya = list(Y)
        #print 'Xa = ', len(Xa), Xa
        #print 'Ya = ', len(Ya), Ya
        #print 'min(Xa) max(Xa) = ', min(Xa), max(Xa)
        #print 'min(Ya) max(Ya) = ', min(Ya), max(Ya)
        #print 'Beam On Off = ', BeamOnOffList
        #print 'min max of BeamOnOffList = ', min(min(BeamOnOffList)), max(max(BeamOnOffList))        
        #print 'New Method ---------------------'
        if min(Xa) > min(min(BeamOnOffList)): 
            Xa.insert(0, min(min(BeamOnOffList)))
            Ya.insert(0, Y[0])
            
        if max(Xa) < max(max(BeamOnOffList)):
            Xa.append(max(max(BeamOnOffList)))
            Ya.append(Y[nX-1])
            
        # print 'Xa = ', len(Xa), Xa
        # print 'Ya = ', len(Ya), Ya
        BState = False 
        nXa = len(Xa)
        nSw = len(BeamOnOffList)
        
        for k in range(nSw):
            XBOnOff, YBOnOff    = [], []       # Empty
            xSw = BeamOnOffList[k]
            for i in range(nXa-1):
                Ysec = Ya[i]
                if Xa[i] <= xSw[0] and xSw[0] <= Xa[i+1]:
                    Ysec = self.lineIntp(xSw[0], [Xa[i], Ya[i]], [Xa[i+1], Ya[i+1]])
                    XBOnOff.append(xSw[0])
                    YBOnOff.append(Ysec)
                    break
            
            for i in range(nXa-1):
                Ysec = Ya[i]
                # print i, Xa[i], xSw[1], Xa[i+1]
                if Xa[i] <= xSw[1] and xSw[1] <= Xa[i+1]:
                    Ysec = self.lineIntp(xSw[1], [Xa[i], Ya[i]], [Xa[i+1], Ya[i+1]])
                    XBOnOff.append(xSw[1])
                    YBOnOff.append(Ysec)
                    break
                    
            #print 'MeanStdEachBeam:: k, XBOnOff YBOnOff = ', k, XBOnOff, YBOnOff
            XBOn, YBOn = [XBOnOff[0]], [YBOnOff[0]]
            for i in range(nX):
                if X[i] > XBOnOff[1]: 
                    break
                #print XBOnOff[0], X[i], XBOnOff[1]
                if XBOnOff[0] < X[i] and X[i] < XBOnOff[1]:
                    XBOn.append(X[i])
                    YBOn.append(Y[i])
                    
            XBOn.append(XBOnOff[1])
            YBOn.append(YBOnOff[1])
            
            if XBOn[0] >= 0.0:
                # print 'MeanStdEachBeam:: k, XBOn YBOn = ', k, XBOn, YBOn
                    
                # print 'self.MeanStdOverall(XBOn, YBOn, M) = ', self.MeanStdOverall(XBOn, YBOn, M)
                #print '[xSw, [YBOn[0], Ysec] = ', xSw, [YBOn[0], Ysec]
                #print 'Ymin Ymax = ', min(YBOn), max(YBOn)
                BeamMeanStd.append([xSw, [YBOn[0], Ysec], self.MeanStdOverall(XBOn, YBOn, M),[min(YBOn), max(YBOn)]])

        # print 'BeamMeanStd = ', BeamMeanStd
        return BeamMeanStd

    def MeanStdSumBeams(self, BeamOnOffMeanStd):
        nBeamOnOffMeanStd = len(BeamOnOffMeanStd)
        if nBeamOnOffMeanStd == 0: return[0,0,0]
        # print "MeanStdSumBeams:: nBeamOnOffMeanStd = ", nBeamOnOffMeanStd
        SumTime = 0
        SumAvg, SumStd = 0, 0
        Rmin, Rmax = [], []
        
        for i in range(nBeamOnOffMeanStd):
            Time = BeamOnOffMeanStd[i][0]  # Time
            Xdat = BeamOnOffMeanStd[i][1]  # Ron and Roff 
            Ydat = BeamOnOffMeanStd[i][2]  # Rmean and Rstddev
            Zdat = BeamOnOffMeanStd[i][3]  # Rmin and Rmax
            Rmin.append(Zdat[0])
            Rmax.append(Zdat[1])
            if Time[0] >= 0.0:
                # print 'MeanStdSumBeams::', i, Time, Xdat, Ydat
                dTime = Time[1] - Time[0]
                # print i, delTime, Ydat[0], Ydat[1]
                SumTime += dTime
                SumAvg  += Ydat[0]*dTime
                # SumStd  += Ydat[1]*Ydat[1]*dTime*dTime/(Ydat[0]*Ydat[0])
                #SumStdX += Std[1]  - Std[0]
                
        m = SumAvg/SumTime
        
        for i in range(nBeamOnOffMeanStd):
            Time = BeamOnOffMeanStd[i][0]  # Time
            Xdat = BeamOnOffMeanStd[i][1]  # Ron and Roff 
            Ydat = BeamOnOffMeanStd[i][2]  # Rmean and Rstddev
            Zdat = BeamOnOffMeanStd[i][3]  # Rmin and Rmax
            Rmin.append(Zdat[0])
            Rmax.append(Zdat[1])
            if Time[0] >= 0.0:
                dTime = Time[1] - Time[0]
                SumStd  += (Ydat[1]*Ydat[1]+(Ydat[0] - m)*(Ydat[0] - m))*dTime
                #SumStdX += Std[1]  - Std[0]       
            
        # print "Sum Time = ", SumTime, " SumAvg = ", SumAvg, " SumStd = ", SumStd
        # print " Avg = ", SumAvg/SumTime, " Std = ", SumStd/SumTime, " Tme = ", SumTime/60.0
        # print '** MeanStdSumBeams:: Rmin, Rmax = ', min(Rmin), max(Rmax)
        # Return Avg, Std, Time in min, Rmin, Rmax
        return [SumAvg/SumTime, NUMPY.sqrt(SumStd/SumTime), SumTime/60.0, min(Rmin), max(Rmax)]
                  
    def updatePlot2dError(self):
        
        iMin, iMax = self.iMin, self.iMax
                        
        T = self.E[iMin:iMax+1]
        R = self.R[iMin:iMax+1]
        X = self.X[iMin:iMax+1]
        Y = self.Y[iMin:iMax+1]
        Z = self.Z[iMin:iMax+1]

        self.plot2dError(X, Y, T, 'XY')  
        self.plot2dError(Y, Z, T, 'YZ')  
        self.plot2dError(X, Z, T, 'XZ')
                
                
    def setEventRange(self):
        
        # print 'itemSelectionChanged ()'
        nRows = self.TableWidgetEvent.rowCount()
        # print 'setEventRange:: nRows = ', nRows
        iMin, iMax = nRows, 0
        for i in range(nRows):
            if self.TableWidgetEvent.item(i,0).isSelected():
                #print i, ' is selected'
                if i <= iMin: iMin = i
                if i >= iMax: iMax = i
                

        if iMin == nRows or iMax == 0: return      
        if iMin == iMax: return          
        print 'setEventRange:: iMin, iMax = ', self.TableWidgetEvent.item(iMin,2).text(), \
                               self.TableWidgetEvent.item(iMax,2).text()
        S = self.TableWidgetEvent.item(iMin,2).text().split(".")
        ms = int(S[1])
        hms = S[0].split(":")
        h = int(hms[0])
        m = int(hms[1])
        s = int(hms[2])
            
        # self.T0 = h*60*60 + m*60 + s + ms/1000.0
        self.T0 = self.time2sec(h,m,s,ms)
        
        
        # Check Reference Points
        #print 'n of Ref Points = ', len(self.RefLines)
        #for line in self.RefLines:
        #    K = list(self.extract_numbers(line))
        #    print 'Ref Point = ', self.time2sec(K[1],K[2],K[3],K[4]), K[7], K[8], K[9]
        
        i = iMin
        self.RefPoint = []
        while (i > 0):
            if "Reference point" in self.TableWidgetEvent.item(i,0).text():
                for line in self.RefLines:
                    K = list(self.extract_numbers(line))
                    if int(self.TableWidgetEvent.item(i,1).text()) == int(K[0]):
                        self.RefPoint = [K[7], K[8], K[9]]
                        self.T0       = self.time2sec(K[1],K[2],K[3],K[4])
                i = 0
            
            i -= 1
        
        # print 'setEventRange:: self.T0 = ', self.T0
        # print 'setEventRange:: self.RefPoint = ', self.RefPoint
        
        self.comboBoxRefTime.clear()
        self.comboBoxRefTime.addItem(str(self.T0))
        
        self.comboBoxRefPoint.clear()
        sItem = ', '.join(str(x) for x in self.RefPoint)
        self.comboBoxRefPoint.addItem(sItem)
        
        sTime = QTime()
        eTime = QTime()
        
        sT = self.TableWidgetEvent.item(iMin,2).text()
        eT = self.TableWidgetEvent.item(iMax,2).text()
        # print 'setEventRange:: s and e = ', sT, eT
        
        S = sT.split(".")
        ms = int(S[1])
        hms = S[0].split(":")
        h = int(hms[0])
        m = int(hms[1])
        s = int(hms[2])
        sTime.setHMS(h, m, s, ms)
        self.timeEditStart.setTime(sTime)
        tMin = h*60*60 + m*60 + s + ms/1000.0
        

        E = eT.split(".")
        ms = int(E[1])
        hms = E[0].split(":")
        h = int(hms[0])
        m = int(hms[1])
        s = int(hms[2])
        eTime.setHMS(h, m, s, ms)
        self.timeEditEnd.setTime(eTime)
        tMax = h*60*60 + m*60 + s + ms/1000.0
        
        if tMin == tMax: return
        
        # Log Number
        lMin = int(self.TableWidgetEvent.item(iMin,1).text())
        lMax = int(self.TableWidgetEvent.item(iMax,1).text())
        # print 'setEventRange:: Log Range = ', lMin, lMax
        self.spinBoxLogMin.setValue(lMin)
        self.spinBoxLogMax.setValue(lMax)
        # Clearance
        # n = len(self.IFMMAlarmLevels)
        # print 'n of IFMM Alarm Levels = ', n
        
        #self.alarmLevel = 1.5
        #for i in range(1,n):
        #    L = self.IFMMAlarmLevels[i-1]
        #    H = self.IFMMAlarmLevels[i]
        #    print i, self.time2sec(L[1],L[2],L[3],L[4]), L[14]
        #    #print i, L[0], H[0]
        #    # if lMin > L[0] and lMin < H[0]: 
        #    if L[0] <= lMax: 
        #        self.alarmLevel = L[14]
        #        #self.alarmLevel = H[14]
        #        #print int(L[0]), ' <= ', lMax
        #        #break
        # 
        #print 'self.alarmLevel = ', self.alarmLevel
            
        # Find iMin and iMax
        iMin, iMax = len(self.L), 0
        for i in range(len(self.L)):
            if self.L[i] <= lMin:
                iMin = i
            else:
                i = len(self.L)
        
        for i in range(iMin, len(self.L)):
            if self.L[i] <= lMax:
                iMax = i
            else:
                i = len(self.L)
                
        # if iMin >= iMax: return    
        self.iMin, self.iMax = iMin, iMax 
        # print 'setEventRange:: Event Range Matrix Index = ', self.iMin, ' to ', self.iMax
        
        T = self.E[iMin:iMax+1]
        R = self.R[iMin:iMax+1]
        X = self.X[iMin:iMax+1]
        Y = self.Y[iMin:iMax+1]
        Z = self.Z[iMin:iMax+1]
        
        # print 'setEventRange:: len of T, R, X, Y, Z = ', len(T), len(R), len(X), len(Y), len(Z)
        if len(T) == 0 or len(R) == 0 or len(X) == 0 or len(Y) == 0 or len(Z) == 0: return
        
        P = NUMPY.sqrt(NUMPY.power(NUMPY.asarray(X)-self.RefPoint[0],2) \
           +NUMPY.power(NUMPY.asarray(Y)-self.RefPoint[1],2) \
           +NUMPY.power(NUMPY.asarray(Z)-self.RefPoint[2],2))
        # print P
        # print 'T = ', T
        # print 'E[iMin] and E[iMax]= ', self.E[iMin], self.E[iMax]
        if self.checkBoxWrite.isChecked():
            HOME  = self.HomeLineEdit.text()
            gkDir = ''
            if self.GKDirListWidget.count() > 0:
                gkDir = self.GKDirListWidget.currentItem().text()
    
            gkLog = self.GKLogListWidget.currentItem().text()
            gkLogFile = gkLog
            if self.checkBoxWriteFile.isChecked():
                gkLogFile = str(DATETIME.date.today())
                            
            RstFile = str(HOME +'/'+ gkDir+'/'+ gkLogFile + '.rst')
            # print 'setEventRange:: Result File = ', RstFile
            txtFileRst = open(RstFile, "a")

            #print RstFile 
            #print gkLog, self.timeEditStart.text(), ' - ', self.timeEditEnd.text(),\
            #      max(X)-min(X), max(Y)-min(Y), max(Z)-min(Z), max(R)-min(R)




            start = str(self.timeEditStart.text())
            end   = str(self.timeEditEnd.text())
            # print 'setEventRange:: start and end = ', start, end
            start_dt = DATETIME.datetime.strptime(start, '%I:%M:%S %p')
            end_dt   = DATETIME.datetime.strptime(end,   '%I:%M:%S %p')
            # print 'setEventRange:: start_dt and end_dt = ', start_dt, end_dt
            diff = (end_dt - start_dt) 
            TxTime = diff.seconds/60.0 
                    
            # print "setEventRange:: max of x, y, z, R = ", max(X), max(Y), max(Z), max(R)
            # print "setEventRange:: min of x, y, z, R = ", min(X), min(Y), min(Z), min(R)
            # print 'RefPoint = ', self.RefPoint
            txtFileRst.write(gkLog+'\t'+self.timeEditStart.text()+' - '+self.timeEditEnd.text() 
                  + '\tREF\t'+str(self.RefPoint[0])+'\t'+str(self.RefPoint[1])+'\t'+str(self.RefPoint[2]) 
                  + '\tTxyzR\t' + str(TxTime) +'\t'
                  + str(max(X)-min(X))+'\t'+str(max(Y)-min(Y))+'\t'+str(max(Z)-min(Z))+'\t'+str(max(R)-min(R)))
            txtFileRst.close()
        #[Rmean, Rstd] = self.MeanStdOverall(T, R)
        #[Xmean, Xstd] = self.MeanStdOverall(T, X)
        #[Ymean, Ystd] = self.MeanStdOverall(T, Y)
        #[Zmean, Zstd] = self.MeanStdOverall(T, Z)
                
        # print 'R: mean +/- std = ', Rmean, Rstd
        # print 'X: mean +/- std = ', Xmean, Xstd
        # print 'Y: mean +/- std = ', Ymean, Ystd
        # print 'Z: mean +/- std = ', Zmean, Zstd
        
        #print 'Avg R while Beam ON = '
        self.AvgStdR = self.MeanStdEachBeam(T, R)
        #print 'Avg X while Beam ON = '
        self.AvgStdX = self.MeanStdEachBeam(T, X)
        #print 'Avg Y while Beam ON = '
        self.AvgStdY = self.MeanStdEachBeam(T, Y)
        #print 'Avg Z while Beam ON = '
        self.AvgStdZ = self.MeanStdEachBeam(T, Z)
        
        self.plot2d(T, R, 'Radial')
        self.plot2dXYZ(T, X, 'X')
        self.plot2dXYZ(T, Y, 'Y')
        self.plot2dXYZ(T, Z, 'Z')
        
        if len(T) < 2: return
        
        self.plot3d(T, X, Y, Z)  
                
            
        self.Xmean, self.Ymean, self.Zmean = 0, 0, 0
        self.Xstd,  self.Ystd,  self.Zstd  = 0, 0, 0
        self.Xdist, self.Ydist, self.Zdist = 0, 0, 0
        
        self.plot2dError(X, Y, T, 'XY')  
        self.plot2dError(Y, Z, T, 'YZ')  
        self.plot2dError(X, Z, T, 'XZ')      
        
        if self.checkBoxWrite.isChecked():
            HOME  = self.HomeLineEdit.text()
            gkDir = ''
            if self.GKDirListWidget.count() > 0:
                gkDir = self.GKDirListWidget.currentItem().text() 
                
            gkLog = self.GKLogListWidget.currentItem().text()
            gkLogFile = gkLog 
            if self.checkBoxWriteFile.isChecked():
                gkLogFile = str(DATETIME.date.today())
                
            RstFile = str(HOME +'/'+ gkDir+'/'+gkLogFile + '.rst')
            txtFileRst = open(RstFile, "a")
            # print Opt, ' mean = ', Xmean, Ymean, '  Distance = ', NUMPY.abs(Xmean-X[0]), NUMPY.abs(Ymean-Y[0])

            if self.checkBoxBeamOverall.isChecked():
                txtFileRst.write('\t'+'XYZ'+'\t'+str(self.Xmean)+'\t'+str(self.Ymean)+'\t'+str(self.Zmean)+'\t'\
                                     +'Std'+'\t'+str(self.Xstd) +'\t'+str(self.Ystd) +'\t'+str(self.Zstd) +'\t'\
                                     +'Dst'+'\t'+str(self.Xdist)+'\t'+str(self.Ydist)+'\t'+str(self.Zdist)+'\n') 
                
            if self.checkBoxBeamOn.isChecked() or self.checkBoxBeamEach.isChecked():
                #nAvgStdR = len(self.AvgStdR)  
                #nAvgStdX = len(self.AvgStdX)  
                #nAvgStdY = len(self.AvgStdY)  
                #nAvgStdZ = len(self.AvgStdZ)          
                #print nAvgStdX, " AvgStdX = ", self.AvgStdX
                # print Opt+" nBeamOnOffMeanStd = ", nBeamOnOffMeanStd
                AvgStdR = self.MeanStdSumBeams(self.AvgStdR)
                AvgStdX = self.MeanStdSumBeams(self.AvgStdX)
                AvgStdY = self.MeanStdSumBeams(self.AvgStdY)
                AvgStdZ = self.MeanStdSumBeams(self.AvgStdZ)
                # print "AvgStdX: XYZ\t", AvgStdX[0],AvgStdY[0],AvgStdZ[0], "\tStd\t", AvgStdX[1],AvgStdY[1],AvgStdZ[1]
                # print "AvgStdR: Rng\t", AvgStdX[3],AvgStdY[3],AvgStdZ[3],AvgStdR[3]
                txtFileRst.write('\t'+'XYZ'+'\t'+str(AvgStdX[0])+'\t'+str(AvgStdY[0])+'\t'+str(AvgStdZ[0]) +'\t'\
                                     +'Std'+'\t'+str(AvgStdX[1])+'\t'+str(AvgStdY[1])+'\t'+str(AvgStdZ[1]) +'\t'\
                                     +'Dst'+'\t'+str(self.Xdist)+'\t'+str(self.Ydist)+'\t'+str(self.Zdist) +'\t'\
                                     +'Tme'+'\t'+str(AvgStdR[2])+'\t' \
                                     +'Rng'+'\t'+str(AvgStdX[4]-AvgStdX[3])+'\t'+str(AvgStdY[4]-AvgStdY[3])+'\t'+str(AvgStdZ[4]-AvgStdZ[3])+'\t'+str(AvgStdR[4]-AvgStdR[3])+'\n') 
 
            txtFileRst.close()
            

        

            

    def getLog(self, logFile, txt):
        
        txtFile = open(logFile, "r")
        
        lineList = []
        logLines = [txt]
        i, k = 0, 0
        
        for line in txtFile:
            lineList.append(line.rstrip("\n"))
            if txt in lineList[i]:    
                K = list(self.extract_numbers(lineList[i]))  
                logLines.append(K)       
                k += 1
                      
            i += 1
        txtFile.close() 
        
        return logLines

                
    def extract_numbers(self, txt):
        for m in re.finditer(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', txt):
            yield float(m.group(0))
        
        
    def readGKLog(self):
        
        # Table Setting
        #if "PyQt5" in sys.modules:
        #    self.TableLabelEvent = QStringListModel()
        #if "PyQt4" in sys.modules:
        #    self.TableLabelEvent = QStringList()
            
        #self.TableLabelEvent << "File Count" << "H" \
        #<< "M"  << "Sec" << "mSec" << "Time(ms)" << "Radial" \
        #<< "x" << "y" << "z" 
        
        self.TableLabelEvent = ["File Count", "H", "M", "Sec", "mSec", "Time(ms)",
                                "Radial", "x", "y", "z"]
            
        home  = self.HomeLineEdit.text()
        gkDir = ''
        if self.GKDirListWidget.count() > 0:
            gkDir = self.GKDirListWidget.currentItem().text() 
            
        gkLog = self.GKLogListWidget.currentItem().text()
        gkLogFile = str(home +'/'+ gkDir + '/' + gkLog + '.log')
#        print '*** gkLogFile = ', gkLogFile
            
        # txtFile = open("C:\PyQt4\ICON\gkcs(2017-01-23 00.01.31).log", "r")
        txtFile = open(gkLogFile, "r")
        
        self.initTableShot()        
              
        lineList = []
        logLine =[]
        i = 0
        k = 0
        sTime = QTime()
        eTime = QTime()

        # List of IFMM Alarm ON  ----------------------
        IMFFAlarmOnList = []          
        for line in txtFile:
            lineList.append(line.rstrip("\n"))
            # New Version of Logfile
            if "IFMM_ALARM_ON_OR_RESTORING" in lineList[i]:
                IMFF = list(self.extract_numbers(lineList[i]))
                IMFFAlarmOnList.append(IMFF)
    
            # Old Version of Logfile
            if "latest marker" in lineList[i]:
                K = list(self.extract_numbers(lineList[i]))  
                logLine.append(K)       
                k += 1
            #          Print result on prompt 
                #print((i+1,k,K))
                T = K[1]*60*60 + K[2]*60 + K[3] + K[4]/1000.0
                # print str(k)+"\t"+str(K[0])+"\t"+str(K[1])+"\t"+str(K[2])+"\t"+str(K[3])+"\t"+str(K[4])+"\t"+str(T)+"\t"+str(K[8])+"\t"+str(K[9])+"\t"+str(K[10])+"\t"+str(K[11])         
            i += 1
        
        #print '*** IMFFAlarmOnList = ', IMFFAlarmOnList
        self.IMFFAlarmOnTime = []
        for K in IMFFAlarmOnList:
            T = K[1]*60*60 + K[2]*60 + K[3] + K[4]/1000.0
            self.IMFFAlarmOnTime.append(T)
            

        # Remove All Rows in the Table
        nRows = self.TableWidgetLog.rowCount()
        for i in range(0, nRows):
            self.TableWidgetLog.removeRow(nRows - i - 1)
        
        # Remove All Columns in the TAble
        nColumns = self.TableWidgetLog.columnCount()
        for i in range(0, nColumns):
            self.TableWidgetLog.removeColumn(nColumns - i - 1)
        
        # Prepare insert rows
        nRows = k
        for i in range(nRows):
            self.TableWidgetLog.insertRow(i)

        #self.setColumnWidth(0,boardSize)
        #self.setColumnWidth(1,85)
        #self.setColumnWidth(2,70)
        #self.setColumnWidth(3,70)
        #self.setColumnWidth(4,70)
        #self.setColumnWidth(5,70)            

        self.TableWidgetLog.setColumnCount(len(self.TableLabelEvent))
        self.TableWidgetLog.setHorizontalHeaderLabels(self.TableLabelEvent) 

        iFile = self.TableLabelEvent.index('File Count')
        iHour = self.TableLabelEvent.index('H')
        iMin  = self.TableLabelEvent.index('M')
        iSec  = self.TableLabelEvent.index('Sec')
        imSec = self.TableLabelEvent.index('mSec')
        iTime = self.TableLabelEvent.index('Time(ms)')
        iRad  = self.TableLabelEvent.index('Radial')
        ix    = self.TableLabelEvent.index('x')
        iy    = self.TableLabelEvent.index('y')
        iz    = self.TableLabelEvent.index('z')        

        k = 0
        self.E, self.R, self.X, self.Y, self.Z = [], [], [], [], []
        self.L = []
        logMin, logMax = 0,0
        # print '*** len(logLine) = ', len(logLine)
        for K in logLine:
            if k == 0:
                # print 'readGKLog:: K[1], K[2], K[3], K[4] = ', K[1], K[2], K[3], K[4]
                sTime.setHMS(K[1], K[2], K[3], K[4])
                # print 'readGKLog:: sTime = ', sTime
                logMin = int(K[0])
                
            T = K[1]*60*60 + K[2]*60 + K[3] + K[4]/1000.0
            
            qTableItem = QTableWidgetItem(str(int(K[0])))
            self.TableWidgetLog.setItem(k,iFile,qTableItem) 
            
            qTableItem = QTableWidgetItem(str(int(K[1])))
            self.TableWidgetLog.setItem(k,iHour,qTableItem) 
            
            qTableItem = QTableWidgetItem(str(int(K[2])))
            self.TableWidgetLog.setItem(k,iMin,qTableItem) 
            
            qTableItem = QTableWidgetItem(str(int(K[3])))
            self.TableWidgetLog.setItem(k,iSec,qTableItem) 
            
            qTableItem = QTableWidgetItem(str(int(K[4])))
            self.TableWidgetLog.setItem(k,imSec,qTableItem) 

            qTableItem = QTableWidgetItem(str(T))
            self.TableWidgetLog.setItem(k,iTime,qTableItem)  
            
            qTableItem = QTableWidgetItem(str(K[8]))
            self.TableWidgetLog.setItem(k,iRad,qTableItem)
            
            qTableItem = QTableWidgetItem(str(K[9]))
            self.TableWidgetLog.setItem(k,ix,qTableItem)
            
            qTableItem = QTableWidgetItem(str(K[10]))
            self.TableWidgetLog.setItem(k,iy,qTableItem)         
            
            qTableItem = QTableWidgetItem(str(K[11]))
            self.TableWidgetLog.setItem(k,iz,qTableItem) 
            
            self.E.append(T)
            self.R.append(K[8])
            self.X.append(K[9])
            self.Y.append(K[10])
            self.Z.append(K[11])
            self.L.append(K[0])  # Log Number
            
            eTime.setHMS(K[1], K[2], K[3], K[4])
            logMax = int(K[0])
            
            k += 1
            # print str(k)+"\t"+str(K[0])+"\t"+str(K[1])+"\t"+str(K[2])+"\t"+str(K[3])+"\t"+str(K[4])+"\t"+str(T)+"\t"+str(K[8])+"\t"+str(K[9])+"\t"+str(K[10])+"\t"+str(K[11])           

        self.timeEditStart.setDisplayFormat('hh:mm:ss AP')
        self.timeEditEnd.setDisplayFormat('hh:mm:ss AP')
        self.timeEditStart.setTime(sTime)
        self.timeEditEnd.setTime(eTime)

        self.spinBoxLogMin.setValue(logMin)
        self.spinBoxLogMax.setValue(logMax)        

        self.TableWidgetLog.resizeColumnsToContents()
        self.TableWidgetLog.resizeRowsToContents()    

        #self.plot2d(self.E, self.R, 'Radial')
        #self.plot2dXYZ(self.E, self.X, 'X')
        #self.plot2dXYZ(self.E, self.Y, 'Y')
        #self.plot2dXYZ(self.E, self.Z, 'Z')
        #self.plot3d(self.X, self.Y, self.Z)
                    

        # print '*** Before readEvents(gkLogFile)'
        self.readEvents(gkLogFile)
        # print '*** After readEvents(gkLogFile)'
        #self.readEvents(gkLogFile, "Powering ECU", "Yellow")
        ## self.readEvents(gkLogFile, "QA_LIST_LOADED")
        #self.readEvents(gkLogFile, "STAND_ALONE_CBCT")
        #self.readEvents(gkLogFile, "EMERGENCY_ALARM_QA")
        #self.readEvents(gkLogFile, "CBCT_PRECISION_QA")
        #self.readEvents(gkLogFile, "FOCUS_PRECISION_QA")
        #self.readEvents(gkLogFile, "TREATMENT_LOADED")
        ## self.readEvents(gkLogFile, "runType: TREATMENT_RUN_DATA_TYPE_TREATMENT")
        ## self.readEvents(gkLogFile, "PERFORMING_TREATMENT")
        #self.readEvents(gkLogFile, "EXECUTING_SHOTS")
        #self.readEvents(gkLogFile, ":CBCT_SCAN")
        #self.readEvents(gkLogFile, "CBCT_DONE")
        ## self.readEvents(gkLogFile, "SENDING_TREATMENT_RUN")
        #self.readEvents(gkLogFile, "SENDING_RUN", 'Green')
        #self.readEvents(gkLogFile, "IDLE")
        ## self.readEvents(gkLogFile, "Leksell DICOM transform")
        #self.readEvents(gkLogFile, "TREATING_PATIENT")
        #self.readEvents(gkLogFile, "EVALUATION_ENDED", 'Red')
        ## self.readEvents(gkLogFile, "ENDING_TREATMENT")
        #self.readEvents(gkLogFile, "PPC1ShotStart")
        #self.readEvents(gkLogFile, "FINALIZING_RUN", 'Blue')
        
        self.plotWidgetXY.canvas.ax.clear()
        self.plotWidgetYZ.canvas.ax.clear()
        self.plotWidgetXZ.canvas.ax.clear()

        
    def readEvents(self, gkLogFile):
        eventList = ["Powering ECU", \
                     "EMERGENCY_ALARM_QA",     "CBCT_PRECISION_QA", \
                     "FOCUS_PRECISION_QA",     "TREATMENT_LOADED", \
                     "SENDING_RUN", "ALL_GOING_TO_HOME_POS",\
                     "CANCELING_RUN",          "ALARM_SET",\
                     # "IN_TRANSPORTATION", "OUT_TRANSPORTATION", \
                     # "ACKNOWLEDGE", \
                     "SCAN_COMPLETED", "SENDING_SCAN_REQUEST",\
                     # "CLOSING", \
                     "TREATMENT_CANCELED_GATING_TIMEOUT", \
                     "IDLE",\
                     "EVALUATION_ACCEPTED",\
                     "treatment timer started", \
                     #"PPC1ShotStart", "STAND_ALONE_CBCT", "CBCT_DONE",\
                     "FINALIZING_RUN", "Reference point", "Sectors beam on state"] 
                
        White  = QColor(255,255,255,255)
        Silver = QColor(169,169,169,127)
        Gray   = QColor(128,128,128,127)
        Black  = QColor(  0,  0,  0,127)
        Red    = QColor(255,  0,  0,127)
        Maroon = QColor(128,  0,  0,127)
        Yellow = QColor(255,255,  0,127)
        Olive  = QColor(128,128,  0,127)
        Lime   = QColor(  0,255,  0, 50)
        Green  = QColor(  0,128,  0,127)
        Aqua   = QColor(  0,255,255,127)
        Teal   = QColor(255,255,  0,127)
        Blue   = QColor(  0,128,128,127)
        Navy   = QColor(  0,  0,128,127)
        Pink   = QColor(255,150,150,127)
        Fuchsia= QColor(255,  0,255,127)
        Purple = QColor(128,  0,128,127) 
        Orange = QColor(255,128,  0,127) 
        
        ### -------- Events Table ----------------------------------
        # Event Table Setting
        #if "PyQt5" in sys.modules:
        #    self.EventTableLabel = QStringListModel()
        #if "PyQt4" in sys.modules:
        #    self.EventTableLabel = QStringList() 
                    
        #self.EventTableLabel << "Event" << "Log" <<  "Time" << "TID"
        self.EventTableLabel = ["Event", "Log",  "Time", "TID"]
        
        self.comboBoxSort.clear()
        self.comboBoxSort.addItems(self.EventTableLabel)
        self.comboBoxSort.setCurrentIndex(1)
                
        # Remove All Rows in the Table
        nRows = self.TableWidgetEvent.rowCount()
        for i in range(0, nRows):
            self.TableWidgetEvent.removeRow(nRows - i - 1)
        
        # Remove All Columns in the TAble
        nColumns = self.TableWidgetEvent.columnCount()
        for i in range(0, nColumns):
            self.TableWidgetEvent.removeColumn(nColumns - i - 1)
        
        self.TableWidgetEvent.setColumnCount(len(self.EventTableLabel))
        self.TableWidgetEvent.setHorizontalHeaderLabels(self.EventTableLabel) 
                
        # Event Table Setting
        #if "PyQt5" in sys.modules:
        #    self.EventTableLabel = QStringListModel()
        #if "PyQt4" in sys.modules:
        #    self.EventTableLabel = QStringList() 
        
        #self.EventTableLabel << "Event" << "Log Number" << "H" \
        #<< "M"  << "Sec" << "mSec" << "Time(ms)" << "TID" 
        #self.EventTableLabel << "Event" << "Log" <<  "Time" << "TID" 
        self.EventTableLabel = ["Event", "Log",  "Time", "TID"] 
        
        iEvent = self.EventTableLabel.index('Event')
        iLog   = self.EventTableLabel.index('Log')
        iTime  = self.EventTableLabel.index('Time')
        iTID   = self.EventTableLabel.index('TID')
       
         
        sTime = QTime()  
        txtFile = open(gkLogFile, "r")
        
        lineList = []
        i, k = 0, 0
        
        self.RefLines = []
        Events        = []
        for line in txtFile:
            lineList.append(line.rstrip("\n"))
            for txt in eventList:
                if txt in lineList[i]:
                    keyword = txt
                    #if txt == "PPC1ShotStart":
                    #    Str     = lineList[i].split(',')
                    #    ShotId  = Str[2].split('=')
                    #    keyword = keyword + ' ('+ShotId[1]+')'
                    if txt == "treatment timer started":
                        Str     = lineList[i].split(',')
                        ShotId  = Str[1].split('=')
                        keyword = '  ShotID ('+ShotId[1]+')'
                        
                    if txt == "Reference point":
                        if not "APPROVED" in lineList[i]:
                            break

                        self.RefLines.append(lineList[i])
                        
                    if txt == "Sectors beam on state":
                        L = list(self.extract_numbers(lineList[i]))
                        S = re.split('\s+', lineList[i])
                        keyword = 'Beam ' + S[12]
                        
                    if txt == "SCAN_COMPLETED":
                        L = list(self.extract_numbers(lineList[i]))
                        S = re.split('\s+', lineList[i])
                        # print S
                        keyCBCT = ''
                        if S[4] == 'mcu::QaCbctPrec':  keyCBCT = ': CBCT QA' 
                        if S[4] == 'mcu::StandAlone':  keyCBCT = ': Stand Alone'
                        if S[4] == 'mcu::DeliveryEv':  keyCBCT = ': Evaluation'
                        keyword = 'CBCT Done' + keyCBCT
                    
                    if txt == "SENDING_SCAN_REQUEST":
                        keyword = 'CBCT Scan'                   
                                        
                        
                    K = [keyword] + list(self.extract_numbers(lineList[i]))
                    Events.append(K)       
                    k += 1
                      
            i += 1
        txtFile.close()
        
        iLine = 0
        V  = lineList[iLine].split(' ')
        if V[0] == '000001': iLine = 4
        V  = lineList[iLine].split(',')
        VC = V[0].split("Software:")
        self.LineEditVersion.setText(VC[1])

        # print '*** Check Point 1, len(Events)', len(Events)        
        for i in range(1,len(Events)):
            # print Event[0], Event[i]       
            K = Events[i]
            iRow = self.TableWidgetEvent.rowCount()
            self.TableWidgetEvent.insertRow(iRow)
            # print 'iRow = ', iRow, Events[i] 
            ColorCode = White
            qTableItem = QTableWidgetItem(K[0])
            if K[0] == "Powering ECU":     ColorCode = Yellow
            if K[0] == "FINALIZING_RUN":   ColorCode = Orange
            if K[0] == "EVALUATION_ACCEPTED": ColorCode = Pink
            if K[0] == "SENDING_RUN":      ColorCode = Green
            if K[0] == "IDLE":             ColorCode = Purple
            if K[0] == "CANCELING_RUN":    ColorCode = Yellow
            if K[0] == "Reference point":  ColorCode = Blue
            if K[0] == "ALARM_SET":        ColorCode = Red
            if K[0] == "TREATMENT_LOADED": ColorCode = Silver
            if K[0] == "ALL_GOING_TO_HOME_POS": ColorCode = Navy
            if 'CBCT' in K[0]:             ColorCode = Lime
            
            qTableItem.setBackground(ColorCode)
            self.TableWidgetEvent.setItem(iRow,iEvent,qTableItem)
            
            # qTableItem.setBackground(White)
            
            qTableItem = QTableWidgetItem(str(int(K[1])))
            self.TableWidgetEvent.setItem(iRow,iLog,qTableItem) 

            sTime.setHMS(K[2], K[3], K[4], K[5])
            qTableItem = QTableWidgetItem(sTime.toString('hh:mm:ss.zzz'))
            self.TableWidgetEvent.setItem(iRow,iTime,qTableItem) 

            qTableItem = QTableWidgetItem(str(int(K[5])))
            self.TableWidgetEvent.setItem(iRow,iTID,qTableItem) 
            
                    
        
        # print '*** Check Point 2'
        nRows = self.TableWidgetEvent.rowCount()
        # print '*** nRows x nColumns = ', nRows, self.TableWidgetEvent.columnCount()
        # if nRows > 0: self.TableWidgetEvent.clearSelection()
        
        #for i in range(nRows):
        #    qTableItem = self.TableWidgetEvent.item(i,iEvent)
        #    #print i, iEvent, qTableItem
        #    self.TableWidgetEvent.setItemSelected(qTableItem, False) # TODO
                      
        # print '*** Check Point 3'
        self.TableWidgetEvent.resizeColumnsToContents()
        self.TableWidgetEvent.resizeRowsToContents()
        # self.TableWidgetEvent.sortItems(2)  
        self.TableWidgetEvent.setColumnWidth(iEvent,200)
        self.TableWidgetEvent.setColumnWidth(iLog,60)
        self.TableWidgetEvent.setColumnWidth(iTID,50)
        # print '*** Check Point 4'

        self.IFMMAlarmLevels      = []
        self.IFMMAlarmTimeList    = []
        self.IFMMAlarmLevelList   = []
        self.SectorStateTimeList  = []
        self.SectorStateValueList = []
        self.BeamStateTimeList    = []
        self.BeamStateValueList   = []
        self.ShotIdStartList      = []
        self.ShotIdStopList       = []
        self.ElapsedTimeList      = []
        self.PlanTimeList         = []
        self.TxTimeStartedList    = []
        self.TxTimeStoppedList    = []
        
        # print '*** Check Point 5'
        nLines = range(len(lineList))
        for i in nLines:
            K = []

            if "ALARM_SET level" in lineList[i]:
                L = list(self.extract_numbers(lineList[i]))
                # print 'L = ', L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0, L[15]
                self.IFMMAlarmTimeList.append(L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0)
                
                if str(self.LineEditVersion.text()).find('11.0') > -1:
                    self.IFMMAlarmLevelList.append(L[16])
                if str(self.LineEditVersion.text()).find('11.1') > -1:
                    self.IFMMAlarmLevelList.append(L[15])
                
            # Not Used --------------------------------------------------------------
            if "SECTOR_STATE_" in lineList[i]:
                L = list(self.extract_numbers(lineList[i]))
                S = re.split('\s+', lineList[i])
                # print 'SECTOR = ', L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0, S[10]
                self.SectorStateTimeList.append(L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0)
                self.SectorStateValueList.append(S[10])   
            # -----------------------------------------------------------------------
                
            if "Sectors beam on state" in lineList[i]:
                L = list(self.extract_numbers(lineList[i]))
                S = re.split('\s+', lineList[i])
                # print 'Beam state = ', L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0, S[12]
                self.BeamStateTimeList.append(L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0)
                self.BeamStateValueList.append(S[12])
                
            if "treatment timer started" in lineList[i]:
                L = list(self.extract_numbers(lineList[i]))
                S = re.split('\s+', lineList[i])
                # print "L = ", L
                # print "S = ", S
                self.PlanTimeList.append(L[7])
                self.TxTimeStartedList.append(L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0)
                # print 'treatment timer started = ', L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0, S[11]
                ShotIdTmp  = S[11].split('=')
                ShotId     = ShotIdTmp[1].split(',')
                # print 'ShotId = ', ShotId[0]
                self.ShotIdStartList.append(ShotId[0])
                
            if "treatment timer stopped" in lineList[i]:
                L = list(self.extract_numbers(lineList[i]))
                S = re.split('\s+', lineList[i])
                # print "L = ", L
                # print "S = ", S
                self.ElapsedTimeList.append(L[7])
                self.TxTimeStoppedList.append(L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0)
                # print 'treatment timer stopped = ', L[1]*60*60 + L[2]*60 + L[3] + L[4]/1000.0, S[11]
                ShotIdTmp  = S[11].split('=')
                ShotId     = ShotIdTmp[1].split(',')
                #self.ShotTimeList.append(S[12])
                # print 'ShotId = ', ShotId[0]
                self.ShotIdStopList.append(ShotId[0])
            
            
                                        
            #if "data {" in lineList[i] or "result {" in lineList[i]:
            #    if "data {" in lineList[i]:   K = list(self.extract_numbers(lineList[i])) 
            #    if "result {" in lineList[i]: K = list(self.extract_numbers(lineList[i-1]))
            #    # print 'K = ', K
            #    IFMMAlarmLevels = []
            #    while (not "checksum" in lineList[i]):
            #        i += 1
            #        # print lineList[i]
            #        if "IFMM ALARM LEVEL" in lineList[i]:
            #            if "event_data" in lineList[i+1]:
            #                L = list(self.extract_numbers(lineList[i+1]))
            #                K.append(L[0])
            #                # print 'K = ', K
            #                IFMMAlarmLevels.append(K)
            #                # print '---', lineList[i+1]                            

            #        if "ifmm_alarm_level_change_data" in lineList[i]:
            #            if "new_alarm_level" in lineList[i+1]:
            #                L = list(self.extract_numbers(lineList[i+1]))
            #                if len(K) > 13:
            #                    K.append(L[0])
            #                    if len(K) > 14:
            #                        K = K[0:14]
            #                        K.append(L[0])
            #                    # print 'L[0] K[14] K = ', L[0], K[14], K
            #                    IFMMAlarmLevels.append(K)
            #                # print '---', lineList[i+1]
                                                        
                            
            #    n = len(IFMMAlarmLevels)
            #    # print 'n = ', n
            #    if n > 0:
            #        self.IFMMAlarmLevels.append(IFMMAlarmLevels[n-1])
            #        # print 'Appended ---', IFMMAlarmLevels[n-1]

        # print 'self.IFMMAlarmLevels = ', self.IFMMAlarmLevels
        # print '*** Check Point 6'
        self.simpleEventView()
        
        nShots = len(self.ShotIdStartList)
        print 'readEvents:: N of Shots = ', nShots
        print 'readEvents:: N of Beam ON and OFF = ', len(self.BeamStateTimeList)
        
        iShotID   = self.TableLabelShot.index('ShotID')
        iPlanned  = self.TableLabelShot.index('Planned')
        iElapsed  = self.TableLabelShot.index('Elapsed')
        iDuration = self.TableLabelShot.index('Duration')
        iStarted  = self.TableLabelShot.index('Started')
        iStopped  = self.TableLabelShot.index('Stopped')
        iBeamOn   = self.TableLabelShot.index('Beam On')
        iBeamOff  = self.TableLabelShot.index('Beam Off')
        iMean     = self.TableLabelShot.index('R_Mean')
        iStdDev   = self.TableLabelShot.index('R_StdDev')
           
           
        for i in range(nShots):
            # print i, self.ShotIdStartList[i], self.TxTimeStartedList[i], self.PlanTimeList[i], \
            # ' -> ', self.ShotIdStopList[i], self.TxTimeStoppedList[i], self.ElapsedTimeList[i]
            
            self.TableWidgetShot.insertRow(i)

            qTableItem = QTableWidgetItem(self.ShotIdStartList[i])
            self.TableWidgetShot.setItem(i,iShotID,qTableItem) 
            
            qTableItem = QTableWidgetItem(str(self.TxTimeStartedList[i]))
            self.TableWidgetShot.setItem(i,iStarted,qTableItem)
            
            qTableItem = QTableWidgetItem(str(self.TxTimeStoppedList[i]))
            self.TableWidgetShot.setItem(i,iStopped,qTableItem)
            
            # qTableItem = QTableWidgetItem(str((self.TxTimeStoppedList[i]-self.TxTimeStartedList[i])*1000))
            # self.TableWidgetShot.setItem(i,iDuration,qTableItem)            
            
            qTableItem = QTableWidgetItem(str(self.PlanTimeList[i]))
            self.TableWidgetShot.setItem(i,iPlanned,qTableItem)
            
            qTableItem = QTableWidgetItem(str(self.ElapsedTimeList[i]))
            self.TableWidgetShot.setItem(i,iElapsed,qTableItem)
                   

                        
        self.TableWidgetShot.resizeColumnsToContents()
        self.TableWidgetShot.resizeRowsToContents()                         
 
    def simpleEventView(self):
        nRows = self.TableWidgetEvent.rowCount()
        iEvent = self.EventTableLabel.index('Event')        
        for i in range(nRows):
            qTableItem = self.TableWidgetEvent.item(i,iEvent)
            if qTableItem:
                # print qTableItem, qTableItem.text()
                if qTableItem.text() == "Beam ON" or qTableItem.text() == "Beam OFF":
                    if self.checkBoxSimple.isChecked():
                        self.TableWidgetEvent.hideRow(i)
                    else:
                        self.TableWidgetEvent.showRow(i)
                        
 
    def readSettingsXml(self):
        self.Home = ""
        self.patient = ""
        self.plan = ""        
        inFile = QFile(self.initFile)
        if not inFile.open(QFile.ReadOnly | QFile.Text):
            if "PyQt5" in sys.modules:
                QMessageBox.warning(self, self.tr("SAX Bookmarks"),
                                    self.tr("Cannot read file {0}:\n{1}.")
                                    .format(self.initFile, inFile.errorString()))
            elif "PyQt4" in sys.modules:
                QMessageBox.warning(self, self.tr("SAX Bookmarks"),
                                    self.tr("Cannot read file %1:\n%2.")
                                    .arg(self.initFile)
                                    .arg(inFile.errorString()))
            else:
                raise SystemError("PyQt4 or PyQt5 need to be imported first") 
                
            return
        
        ok, errorStr, errorLine, errorColumn = self.domDocument.setContent(inFile, True)
        
        if not ok:
            if "PyQt5" in sys.modules:
                QMessageBox.information(self.window(), self.tr("DOM Bookmarks"),
                            self.tr("Parse error at line {0}, column {1}:\n{2}")
                            .format(errorLine, errorColumn, errorStr))
            elif "PyQt4" in sys.modules:
                QMessageBox.information(self.window(), self.tr("DOM Bookmarks"),
                              self.tr("Parse error at line %1, column %2:\n%3")
                              .arg(errorLine).arg(errorColumn).arg(errorStr))  
            else:
                raise SystemError("PyQt4 or PyQt5 need to be imported first")

            return False

        inFile.close()
        
        root = self.domDocument.documentElement()
        
        child = root.firstChild()
        while not child.isNull():
            element = child.toElement()
            if element.isElement():
                # print element.tagName()
                if (element.tagName() == "HOME"):
                    self.Home = str(element.attribute("Path", "")).strip()
                    
                if (element.tagName() == "Recent"):
                    self.patient = str(element.attribute("Directory", "")).strip()
                    self.plan = str(element.attribute("Log", "")).strip()
                    
                
            child = child.nextSibling()
            
        if self.Home is not "":
            self.HomeLineEdit.setText(self.Home)

             
        #if self.patient is not "":
        self.updateGKDirList()
            
        #if self.plan is not "":
        self.updateGKLogList()       
        
 
    def saveSettingsXml(self):
        self.domDocument.clear()
        
        self.PlanViewerRoot = self.domDocument.createElement("IconGating")
        self.domDocument.appendChild(self.PlanViewerRoot)

        self.HomeElement = self.domDocument.createElement("HOME")
        self.PlanViewerRoot.appendChild(self.HomeElement)
        self.HomeElement.setAttribute("Path", self.HomeLineEdit.text())

        self.RecentElement = self.domDocument.createElement("Recent")
        self.PlanViewerRoot.appendChild(self.RecentElement)
        self.RecentElement.setAttribute("Directory", self.GKDirListWidget.currentItem().text())  
        self.RecentElement.setAttribute("Log", self.GKLogListWidget.currentItem().text())
         
             
        outFile = QFile(self.initFile)
        if not outFile.open(QFile.WriteOnly | QFile.Text):
            if "PyQt5" in sys.modules:
                QMessageBox.warning(self, self.tr("SAX Bookmarks"),
                                    self.tr("Cannot write file {0}:\n{1}.")
                                    .format(self.initFile, outFile.errorString()))
            elif "PyQt4" in sys.modules:
                QMessageBox.warning(self, self.tr("SAX Bookmarks"),
                                    self.tr("Cannot write file %1:\n%2.")
                                    .arg(self.initFile)
                                    .arg(outFile.errorString()))
            else:
                raise SystemError("PyQt4 or PyQt5 need to be imported first")
                      
            return        
        out = QTextStream(outFile)
        self.domDocument.save(out, 4)
        
        outFile.close()
       

    def setDatabaseHome(self):
        xHome = self.DatabaseHomeLineEdit.text()
        if xHome is '':
            xHome = self.ROOT
        self.DatabaseHome = QFileDialog.getExistingDirectory(self, "Open Database Home...",
            xHome, QFileDialog.ShowDirsOnly)
        if self.DatabaseHome is "" :
            self.DatabaseHome = xHome
        self.DatabaseHomeLineEdit.setText(self.DatabaseHome)
        
        

    def setHome(self):
        xHome = self.HomeLineEdit.text()
        if xHome is '':
            xHome = self.ROOT
        self.Home = QFileDialog.getExistingDirectory(self, "Open Leksell ICON Log File Home...",
            xHome, QFileDialog.ShowDirsOnly)
        if self.Home is "" :
            self.Home = xHome
        self.HomeLineEdit.setText(self.Home)
        
        self.updateGKDirList() 
        
        

        
    def updateGKDirList(self):
        if self.Home == "":
            print "No self.Home assigned..."
            return
            
        os.chdir(str(self.Home))
        CWD = os.getcwd()
        if os.path.exists(CWD):
            os.chdir(CWD)
            gkcsDirList = os.listdir(CWD)
            # self.messageLabel.setText("Searching Patient Directories...")

            self.GKDirListWidget.clear()
            
            gkcsList = []
            for gkcsDir in gkcsDirList:
                if os.path.isdir(gkcsDir):
                    gkcs = os.path.basename(gkcsDir)
                    # print 'gkcs Dir = ', gkcs
                    if gkcsList.count(gkcs) == 0:
                        gkcsList.append(gkcs)
                        
            gkcsList.sort()
            # print patientList
            self.GKDirListWidget.addItems(gkcsList)
            
            # Set Cursor to Previous Patient
            if self.patient is "":
                return
            
            nRows = self.GKDirListWidget.count()
            iCurrentRow = 0
            for i in range(nRows):
                if self.patient == str(self.GKDirListWidget.item(i).text()).strip():
                    iCurrentRow = i
            
            self.GKDirListWidget.setCurrentRow(iCurrentRow) 
           
        else:
            print "No " + CWD
            return
          
        self.updateGKLogList()
        
        
    def updateGKLogList(self):
        
        subDir = ''
        if self.GKDirListWidget.count() > 0:
            subDir = self.GKDirListWidget.currentItem().text()
            
        gkDir = str(self.HomeLineEdit.text()).strip() + "/" + subDir
        # print 'gkDir = ', gkDir
        
        os.chdir(str(gkDir))
        CWD = os.getcwd()
        if os.path.exists(CWD):
            os.chdir(CWD)
            gkLogDirList = os.listdir(CWD)
            # self.messageLabel.setText("Searching Patient Directories...")
            self.GKLogListWidget.clear() 
            
            # -----------------------------------------
            gkLogList = []
            for gkLogDir in gkLogDirList:
                if os.path.isfile(gkLogDir):
                    gkLog = ''
                    gkLogFile = os.path.basename(gkLogDir)
                    if gkLogFile.find(".log") > -1:
                        gkLog = gkLogFile.split(".log")[0]
                        
                        if gkLogList.count(gkLog) == 0:
                            gkLogList.append(gkLog)
                                                     
                if os.path.isdir(gkLogDir):
                    gkLogSubDirList = os.listdir(CWD+'/'+gkLogDir)
                    for gkLogSubDir in gkLogSubDirList:
                        if os.path.isfile(gkLogDir+'/'+gkLogSubDir):
                            gkLog = ''
                            gkLogFile = os.path.basename(gkLogSubDir)
                            if gkLogFile.find(".log") > -1:
                                gkLog = gkLogFile.split(".log")[0]
                                
                                if gkLogList.count(gkLog) == 0:
                                    gkLogList.append(gkLogDir+'/'+gkLog)
            

            gkLogList.sort()
            # print 'gkLogList = ', gkLogList
            self.GKLogListWidget.addItems(gkLogList)
 
            # Set Cursor to Previous Plan
            nItems = self.GKLogListWidget.count()
            if self.plan is "":
                return
            
            iCurrentItem = 0
            for i in range(nItems):
                if str(self.plan).strip() == str(self.GKLogListWidget.item(i).text()).strip():
                    iCurrentItem = i
            
            self.GKLogListWidget.setCurrentRow(iCurrentItem) 
            

                    
    def plot2d(self, X, Y, Label):   
         
        if X == []: 
            self.plotWidget.canvas.ax.clear()
            return     
        
        # Elapse Time from Reference Time
        # XR = NUMPY.asarray(X)-X[0]
        XR = NUMPY.asarray(X)-self.T0
        #print '* plot2d:: XR = ', XR 
        #print '* plot2d:: len(XR) = ', len(XR)
        
        BeamOnOffMeanStd = self.MeanStdEachBeam(X, Y)
        nBeamOnOffMeanStd = len(BeamOnOffMeanStd)
        #print nBeamOnOffMeanStd, " R plot2d::BeamOnOffMeanStd = ", BeamOnOffMeanStd
        #print 'R nBeamOnOffMeanStd = ', nBeamOnOffMeanStd
         
        MeanStdSumBeams = self.MeanStdSumBeams(BeamOnOffMeanStd)
        # print "R Avg = ", MeanStdSumBeams[0], " Std = ", MeanStdSumBeams[1]

                
        self.plotWidget.canvas.ax.clear()
        
        # IMFFAlarmTime = NUMPY.asarray(self.IMFFAlarmOnTime)-X[0]
        IMFFAlarmTime = NUMPY.asarray(self.IMFFAlarmOnTime)-self.T0
        # print 'IMFFAlarmTime = ', IMFFAlarmTime
        for t in IMFFAlarmTime:
            if t > XR[0] and t < max(XR) and self.checkBoxPaused.isChecked():
                self.plotWidget.canvas.ax.axvline(x=t, linewidth=1, color = 'yellow')
                        
        self.plotWidget.canvas.ax.plot(XR, Y, 'o-', color = 'blue', ms=2)
        if Label == 'Radial':
            # self.plotWidget.canvas.ax.axhline(y=self.alarmLevel, linewidth=1, color = 'red')
            nLevel = len(self.IFMMAlarmTimeList)
            # self.IFMMAlarmTimeList.append(0.0)
            # self.IFMMAlarmLevelList.append(1.5)
            IFMMAlarmTimeList = NUMPY.asarray(self.IFMMAlarmTimeList)-self.T0
            IFMMAlarmLevelList = NUMPY.asarray(self.IFMMAlarmLevelList)
            
            xLim = self.plotWidget.canvas.ax.get_xlim()
            IFMMAlarmTimeList  = NUMPY.append(IFMMAlarmTimeList, xLim[1])
            IFMMAlarmLevelList = NUMPY.append(IFMMAlarmLevelList, 1.5)

            # print IFMMAlarmLevelList
            for i in range(nLevel+1):
                print 'plot2d:: Alarm List = ', i, IFMMAlarmTimeList[i], IFMMAlarmLevelList[i]            
                
            self.plotWidget.canvas.ax.step(IFMMAlarmTimeList, IFMMAlarmLevelList, linewidth=1, color = 'red', where = 'post')
            
            nSectorValue = len(self.SectorStateValueList)
            SectorStateTimeList = NUMPY.asarray(self.SectorStateTimeList)-self.T0
            for i in range(nSectorValue):
                #if self.SectorStateValueList[i].find('ALL_AT_TREAT_POS') > 1:
                #    if (SectorStateTimeList[i] > XR[0] and SectorStateTimeList[i] <= max(XR)):
                #        print SectorStateTimeList[i], self.SectorStateValueList[i]
                #        #self.plotWidget.canvas.ax.axvline(x=SectorStateTimeList[i], linewidth=1, color = 'blue')

                #if self.SectorStateValueList[i].find('ALL_GOING_TO_BEAM_OFF_POS') > 1:
                #    if (SectorStateTimeList[i] > XR[0] and SectorStateTimeList[i] <= max(XR)):
                #        print SectorStateTimeList[i], self.SectorStateValueList[i]
                #        #self.plotWidget.canvas.ax.axvline(x=SectorStateTimeList[i], linewidth=1, color = 'red') 
                        
                if self.SectorStateValueList[i].find('ALL_GOING_TO_HOME_POS') > 1:
                    if (SectorStateTimeList[i] > XR[0] and SectorStateTimeList[i] <= max(XR)):
                        # print SectorStateTimeList[i], self.SectorStateValueList[i]
                        self.plotWidget.canvas.ax.axvline(x=SectorStateTimeList[i], linewidth=1, color = 'green')    
            
            # nBeamState = len(self.BeamStateValueList)
            # print 'nBeamState = ', nBeamState
            # BeamStateTimeList = NUMPY.asarray(self.BeamStateTimeList)-self.T0
            # BeamOnOffList = []
            # BeamOn, BeamOff = 0.0, 0.0
            # for i in range(nBeamState):
            #     if (BeamStateTimeList[i] > XR[0] and BeamStateTimeList[i] <= max(XR)):
            #         # print BeamStateTimeList[i], self.BeamStateValueList[i], self.BeamStateValueList[i].find('ON')
            #         if self.BeamStateValueList[i].find('ON')  == 0:
            #             self.plotWidget.canvas.ax.axvline(x=BeamStateTimeList[i], linewidth=1, color = 'blue')
            #             BeamOn = BeamStateTimeList[i]
            #             
            #         if self.BeamStateValueList[i].find('OFF') == 0:
            #             self.plotWidget.canvas.ax.axvline(x=BeamStateTimeList[i], linewidth=1, color = 'red')
            #             BeamOff = BeamStateTimeList[i]
            #             if BeamOff > BeamOn: 
            #                 BeamOnOffList.append([BeamOn, BeamOff])
            # 
            # for x in BeamOnOffList:
            #     print 'Beam On to Off = ', x
            #     iOn, iOff = 0, 0 
            #     XBOnOff, YBOnOff = [], []
            #     nXR = len(XR)-1
            #     for k in range(nXR):
            #         if XR[k] <= x[0] and x[0] < XR[k+1]:
            #             iOn = k
            #             Ysec = self.lineIntp(x[0], [XR[k], Y[k]], [XR[k+1], Y[k+1]])
            #             print 'Beam ON  [',XR[k], Y[k],']', ' <= [', x[0],Ysec, '] < ', '[',XR[k+1], Y[k+1],']'
            #             XBOnOff.append(x[0])
            #             YBOnOff.append(Ysec)
            #             
            #         if XR[k] <= x[1] and x[1] < XR[k+1]:
            #             iOff = k+1
            #             Ysec = self.lineIntp(x[1], [XR[k], Y[k]], [XR[k+1], Y[k+1]])
            #             print 'Beam OFF [',XR[k], Y[k],']', ' <= [', x[1],Ysec, '] < ', '[',XR[k+1], Y[k+1],']'
            #             XBOnOff.append(x[1])
            #             YBOnOff.append(Ysec)
            #             
            #     print 'Beam On to Off Index = ', iOn, iOff+1
            #     Xb = XR[iOn+1:iOff].tolist()
            #     Xb.insert(0, XBOnOff[0])
            #     Xb.append(XBOnOff[1])
            #     
            #     Yb = Y[iOn+1:iOff]
            #     Yb.insert(0, YBOnOff[0])
            #     Yb.append(YBOnOff[1])
            #     
            #     print 'XBOnOff = ', XBOnOff, " YBOnOff = ", YBOnOff, " Mean, Std = ", self.MeanStdOverall(Xb, Yb) 
            #     # print 'X = ', XB[0], XR[iOn:iOff+1], XB[1]
            #     # print 'Y = ', YB[0],  Y[iOn:iOff+1], YB[1]
            #  
            #     # print 'Xb = ', Xb
            #     # print 'Yb = ', Yb
        
        TxTimeStartList   = NUMPY.asarray(self.TxTimeStartedList) - self.T0
        TxTimeStoppedList = NUMPY.asarray(self.TxTimeStoppedList) - self.T0
        # print 'nBeamStartList = ', len(TxTimeStartList)
        # print 'nBeamStoppedList = ', len(TxTimeStoppedList)
        
        iMean     = self.TableLabelShot.index('R_Mean')
        iStdDev   = self.TableLabelShot.index('R_StdDev') 
        iBeamOn   = self.TableLabelShot.index('Beam On')
        iBeamOff  = self.TableLabelShot.index('Beam Off')
        iDuration = self.TableLabelShot.index('Duration')
        
        # Remove All Rows in the Table
        # nRows = self.TableWidgetShot.rowCount()
        # for i in range(0, nRows):
        #     self.TableWidgetShot.removeCellWidget(i,iBeamOn)
        #     self.TableWidgetShot.removeCellWidget(i,iBeamOff)
        #     self.TableWidgetShot.removeCellWidget(i,iDuration)
        #     self.TableWidgetShot.removeCellWidget(i,iMean)
        #     self.TableWidgetShot.removeCellWidget(i,iStdDev)

                    
        nTxTime = len(TxTimeStartList)
        for iBeamOnOff in BeamOnOffMeanStd:
            xOnOff  = iBeamOnOff[0]
            yOnOff  = iBeamOnOff[1]
            meanStd = iBeamOnOff[2]
            
            nList = len(self.TxTimeStartedList)
            for i in range(nList):
                # print 'xOnOff[0]+T0 = ', float(xOnOff[0])+self.T0, self.TxTimeStartedList[i], float(xOnOff[1])+self.T0, self.TxTimeStoppedList[i]
                if float(xOnOff[0])+self.T0 == self.TxTimeStartedList[i] \
                   or float(xOnOff[1])+self.T0 == self.TxTimeStoppedList[i]:
                    
                    # print i, self.TxTimeStartedList[i], self.TxTimeStoppedList[i], meanStd
                    
                    qTableItem = QTableWidgetItem(str(xOnOff[0]))
                    self.TableWidgetShot.setItem(i,iBeamOn,qTableItem)
                    
                    qTableItem = QTableWidgetItem(str(xOnOff[1]))
                    self.TableWidgetShot.setItem(i,iBeamOff,qTableItem)
                    
                    qTableItem = QTableWidgetItem(str((xOnOff[1] - xOnOff[0])/60.0))
                    self.TableWidgetShot.setItem(i,iDuration,qTableItem) 
            
                    qTableItem = QTableWidgetItem(str(meanStd[0]))
                    self.TableWidgetShot.setItem(i,iMean,qTableItem)
                    
                    qTableItem = QTableWidgetItem(str(meanStd[1]))
                    self.TableWidgetShot.setItem(i,iStdDev,qTableItem)
            
                    break
            
            if self.checkBoxBeamEach.isChecked():
                if len(xOnOff) > 2: xOnOff = xOnOff[0:2] # Added
                if len(yOnOff) > 2: yOnOff = yOnOff[0:2] # Added    
                # print 'xOnOff = ', xOnOff, ' yOnOff = ', yOnOff, ' Mean and Std = ', meanStd
                self.plotWidget.canvas.ax.plot(xOnOff, [meanStd[0],meanStd[0]], '-', color = 'm', ms=1)

            
            if self.checkBoxBeamOnOff.isChecked():
                for i in range(nTxTime):
                    #if xOnOff[0] >= TxTimeStartList[i] and xOnOff[0] < TxTimeStoppedList[i]:
                    #    # print 'Shot ID = ', self.ShotIdStartList[i]
                    #    break
                    self.plotWidget.canvas.ax.axvline(x=TxTimeStartList[i],   linewidth=0.5, color = 'blue')
                    self.plotWidget.canvas.ax.axvline(x=TxTimeStoppedList[i], linewidth=0.5, color = 'red')

        if self.checkBoxBeamOn.isChecked():
            BeamOne = BeamOnOffMeanStd[0]
            TimeOne = BeamOne[0]
            for BeamOne in BeamOnOffMeanStd:
                TimeOne = BeamOne[0]
                if TimeOne[0] >= 0.0:
                    break

            BeamEnd = BeamOnOffMeanStd[nBeamOnOffMeanStd-1]            
            TimeEnd = BeamEnd[0]
            # print "TimeOne and TimeEnd = ", TimeOne[0], TimeEnd[1]
            # print "self.checkBoxBeamOn.isChecked(): ", TimeOne[0], TimeEnd[1], MeanStdSumBeams[0], MeanStdSumBeams[1]
            self.plotWidget.canvas.ax.plot([TimeOne[0], TimeEnd[1]], [MeanStdSumBeams[0], MeanStdSumBeams[0]], '-', color = 'm', ms=1)
                
        labels_x = self.plotWidget.canvas.ax.get_xticklabels()
        labels_y = self.plotWidget.canvas.ax.get_yticklabels()
        for xlabel in labels_x:
            xlabel.set_fontsize(7)

        for ylabel in labels_y:
            ylabel.set_fontsize(7)
                        
        xLim = self.plotWidget.canvas.ax.get_xlim()
        yLim = self.plotWidget.canvas.ax.get_ylim()
        # print xLim, yLim
        # zFactor = self.ZoomScrollBar.value()
        #print 'xMin, xMax = ', min(XR), max(XR)
        #print 'x ranges = ', min(XR) - min(XR) % 60, max(XR) - max(XR) % 60
        
        zFactor = 100
        xMin, xMax = 100.0/zFactor*xLim[0], 100.0/zFactor*xLim[1]
        yMin, yMax = 100.0/zFactor*yLim[0], 100.0/zFactor*yLim[1]
        
        xMin, xMax = min(XR), max(XR)
        # print 'max(XR) % 60, max(XR) / 60 = ', min(XR) - min(XR) % 60, max(XR) - max(XR) % 60 + 60.1
        xRanges = NUMPY.arange(min(XR) - min(XR) % 60, max(XR) - max(XR) % 60 + 60, 60)
        # print 'xRanges = ', xRanges

        if xMin < 0: 
            xMin = 0.0
            xRanges = NUMPY.arange(xMin, max(XR) - max(XR) % 60 + 60, 60)

        if self.checkBoxAxisSetting.isChecked():
            xMin = self.spinBoxMinX.value() 
            xMax = self.spinBoxMaxX.value()
            # xMax = xMin + 600
            xRanges = NUMPY.arange(xMin + 60 - xMin % 60, xMax + 60 - xMax % 60, 60)
        
        
        self.spinBoxMinX.setValue(xMin)
        self.spinBoxMaxX.setValue(xMax)
                    
        yMin, yMax = 0.0, max(self.IFMMAlarmLevelList)+0.5
        self.plotWidget.canvas.ax.set_xlim(xMin,xMax)
        self.plotWidget.canvas.ax.set_ylim(yMin,yMax)
        # self.plotWidget.canvas.ax.xaxis.set_ticks(NUMPY.arange(xMin, xMax, 60))
        self.plotWidget.canvas.ax.xaxis.set_ticks(xRanges)
        

        self.plotWidget.canvas.ax.xaxis.set_label_coords(0.5, -0.12) 
        self.plotWidget.canvas.ax.set_xlabel("T(s)", fontsize = 8) 
        self.plotWidget.canvas.ax.set_ylabel("RMS(mm)", fontsize = 8)         
        
        self.plotWidget.canvas.ax.grid(True)
        
        self.plotWidget.canvas.draw()
        
        
    def plot2dXYZ(self, X, Y, Opt='X'):    
            
        if X == []: 
            self.plotWidgetX.canvas.ax.clear()
            self.plotWidgetY.canvas.ax.clear()
            self.plotWidgetZ.canvas.ax.clear()
            return     
                
        # XR = NUMPY.asarray(X)-X[0]
        XR = NUMPY.asarray(X)-self.T0
        plotWidget = self.plotWidgetX

        Yo = self.RefPoint[0]
                
        yLabel = 'X(mm)'
        meanLabel   = Opt+'_Mean'
        stddevLabel = Opt+'_StdDev'
        if Opt == 'X':
            plotWidget = self.plotWidgetX
            yLabel = 'X(mm)'
            Yo = self.RefPoint[0]
            
        if Opt == 'Y':
            plotWidget = self.plotWidgetY
            yLabel = 'Y(mm)'
            Yo = self.RefPoint[1]
            
        if Opt == 'Z':
            plotWidget = self.plotWidgetZ
            yLabel = 'Z(mm)'
            Yo = self.RefPoint[2]
            

        BeamOnOffMeanStd = self.MeanStdEachBeam(X, Y)
        nBeamOnOffMeanStd = len(BeamOnOffMeanStd)            
        # print nBeamOnOffMeanStd, " "+Opt+" BeamOnOffMeanStd = ", BeamOnOffMeanStd
        # print Opt+" nBeamOnOffMeanStd = ", nBeamOnOffMeanStd
        
        MeanStdSumBeams = self.MeanStdSumBeams(BeamOnOffMeanStd)
        # print Opt + " Avg = ", MeanStdSumBeams[0], " Std = ", MeanStdSumBeams[1]
                
        plotWidget.canvas.ax.clear()   
        plotWidget.canvas.ax.plot(XR, Y, 'o-', color = 'blue', ms=2)
        
        # print 'self.BeamStateTimeList = ', self.BeamStateTimeList
        # print 'self.BeamStateValueList = ', self.BeamStateValueList
        # print 'plot2dXYZ:: NUMPY.float64(Yo) = ', NUMPY.float64(Yo)
        [Ym, Ys] = self.MeanStdOverall(X, Y, NUMPY.float64(Yo))
        plotWidget.canvas.ax.axhline(y=Yo, linewidth=1, color = 'green')  
        #plotWidget.canvas.ax.axhline(y=Yo+Ys, linewidth=1, linestyle='--', color = 'red') 
        #plotWidget.canvas.ax.axhline(y=Yo-Ys, linewidth=1, linestyle='--', color = 'red') 
        
        if self.checkBoxBeamOverall.isChecked():
                     
            [Ym, Ys] = self.MeanStdOverall(X, Y)
            plotWidget.canvas.ax.axhline(y=Ym, linewidth=1, color = 'red')     
            #plotWidget.canvas.ax.axhline(y=Ym+Ys, linewidth=1, linestyle='--', color = 'red') 
            #plotWidget.canvas.ax.axhline(y=Ym-Ys, linewidth=1, linestyle='--', color = 'red')        

        # ---------------------------------------------------------------
        TxTimeStartList   = NUMPY.asarray(self.TxTimeStartedList) - self.T0
        TxTimeStoppedList = NUMPY.asarray(self.TxTimeStoppedList) - self.T0
        # print 'nBeamStartList = ', len(TxTimeStartList)
        # print 'nBeamStoppedList = ', len(TxTimeStoppedList)
        
        iMean     = self.TableLabelShot.index(meanLabel)
        iStdDev   = self.TableLabelShot.index(stddevLabel) 
        iBeamOn   = self.TableLabelShot.index('Beam On')
        iBeamOff  = self.TableLabelShot.index('Beam Off')
        iDuration = self.TableLabelShot.index('Duration')
        
        # Remove All Rows in the Table
        nRows = self.TableWidgetShot.rowCount()
        for i in range(0, nRows):
            self.TableWidgetShot.removeCellWidget(i,iBeamOn)
            self.TableWidgetShot.removeCellWidget(i,iBeamOff)
            self.TableWidgetShot.removeCellWidget(i,iDuration)
            self.TableWidgetShot.removeCellWidget(i,iMean)
            self.TableWidgetShot.removeCellWidget(i,iStdDev)
            
                    
        nTxTime = len(TxTimeStartList)
        for iBeamOnOff in BeamOnOffMeanStd:
            xOnOff  = iBeamOnOff[0]
            yOnOff  = iBeamOnOff[1]
            meanStd = iBeamOnOff[2]
            
            nList = len(self.TxTimeStartedList)
            for i in range(nList):
                # print 'xOnOff[0]+T0 = ', float(xOnOff[0])+self.T0, self.TxTimeStartedList[i], float(xOnOff[1])+self.T0, self.TxTimeStoppedList[i]
                if float(xOnOff[0])+self.T0 == self.TxTimeStartedList[i] \
                   or float(xOnOff[1])+self.T0 == self.TxTimeStoppedList[i]:
                    
                    # print i, self.TxTimeStartedList[i], self.TxTimeStoppedList[i], meanStd
                    
                    qTableItem = QTableWidgetItem(str(xOnOff[0]))
                    self.TableWidgetShot.setItem(i,iBeamOn,qTableItem)
                    
                    qTableItem = QTableWidgetItem(str(xOnOff[1]))
                    self.TableWidgetShot.setItem(i,iBeamOff,qTableItem)
                    
                    qTableItem = QTableWidgetItem(str((xOnOff[1] - xOnOff[0])/60.0))
                    self.TableWidgetShot.setItem(i,iDuration,qTableItem) 
            
                    # X_Mean Y_Mean Z_Mean
                    # qTableItem = QTableWidgetItem(str(meanStd[0]))
                    qTableItem = QTableWidgetItem(str(Yo))
                    self.TableWidgetShot.setItem(i,iMean,qTableItem)
                    
                    # X_StdDev Y_StdDev Z_StdDev
                    #qTableItem = QTableWidgetItem(str(meanStd[1]))
                    qTableItem = QTableWidgetItem(str(meanStd[0]-Yo))
                    self.TableWidgetShot.setItem(i,iStdDev,qTableItem)
            
                    break
            
            if self.checkBoxBeamEach.isChecked():
                # print 'xOnOff = ', xOnOff, ' yOnOff = ', yOnOff, ' Mean and Std = ', meanStd
                if len(xOnOff) > 2: xOnOff = xOnOff[0:2] # Added
                if len(yOnOff) > 2: yOnOff = yOnOff[0:2] # Added
                plotWidget.canvas.ax.plot(xOnOff, [meanStd[0],meanStd[0]], '-', color = 'm', ms=1)
            
            if self.checkBoxBeamOnOff.isChecked():
                for i in range(nTxTime):
                    #if xOnOff[0] >= TxTimeStartList[i] and xOnOff[0] < TxTimeStoppedList[i]:
                    #    # print 'Shot ID = ', self.ShotIdStartList[i]
                    #    break
                    
                    plotWidget.canvas.ax.axvline(x=TxTimeStartList[i],   linewidth=0.5, color = 'blue')
                    plotWidget.canvas.ax.axvline(x=TxTimeStoppedList[i], linewidth=0.5, color = 'red')

            
        if self.checkBoxBeamOn.isChecked():
            BeamOne = BeamOnOffMeanStd[0]
            TimeOne = BeamOne[0]
            for BeamOne in BeamOnOffMeanStd:
                TimeOne = BeamOne[0]
                if TimeOne[0] >= 0.0:
                    break

            BeamEnd = BeamOnOffMeanStd[nBeamOnOffMeanStd-1]            
            TimeEnd = BeamEnd[0]
                        
            # print "checkBoxBeamOn.isChecked(): ", TimeOne[0], TimeEnd[1], MeanStdSumBeams[0], MeanStdSumBeams[1]
            plotWidget.canvas.ax.plot([TimeOne[0], TimeEnd[1]], [MeanStdSumBeams[0], MeanStdSumBeams[0]], '-', color = 'm', ms=1)
        # ---------------------------------------------------------------
            
        labels_x = plotWidget.canvas.ax.get_xticklabels()
        labels_y = plotWidget.canvas.ax.get_yticklabels()
        for xlabel in labels_x:
            xlabel.set_fontsize(7)

        for ylabel in labels_y:
            ylabel.set_fontsize(7)
                        
        xLim = plotWidget.canvas.ax.get_xlim()
        yLim = plotWidget.canvas.ax.get_ylim()
        # print xLim, yLim
        # zFactor = self.ZoomScrollBar.value()

        zFactor = 100
        xMin, xMax = 100.0/zFactor*xLim[0], 100.0/zFactor*xLim[1]
        yMin, yMax = 100.0/zFactor*yLim[0], 100.0/zFactor*yLim[1]
        
        xMin, xMax = min(XR), max(XR)
        # print 'max(XR) % 60, max(XR) / 60 = ', min(XR) - min(XR) % 60, max(XR) - max(XR) % 60 + 60.1
        xRanges = NUMPY.arange(min(XR) - min(XR) % 60, max(XR) - max(XR) % 60 + 60, 60)
        # print 'xRanges = ', xRanges
        
        if xMin < 0: 
            xMin = 0.0
            xRanges = NUMPY.arange(xMin, max(XR) - max(XR) % 60 + 60, 60)
                     
        if self.checkBoxAxisSetting.isChecked():
            xMin = self.spinBoxMinX.value() 
            # xMax = self.spinBoxMaxX.value()
            xMax = xMin + 600
            xRanges = NUMPY.arange(xMin + 60 - xMin % 60, xMax + 60 - xMax % 60, 60)
        
        plotWidget.canvas.ax.set_xlim(xMin,xMax)
        plotWidget.canvas.ax.set_ylim(yMin,yMax)
        plotWidget.canvas.ax.xaxis.set_ticks(xRanges)
                
        plotWidget.canvas.ax.xaxis.grid(True)    
         
        plotWidget.canvas.ax.xaxis.set_label_coords(0.5, -0.12) 
        plotWidget.canvas.ax.set_xlabel("T(s)", fontsize = 8) 
        plotWidget.canvas.ax.set_ylabel(yLabel,  fontsize = 8)
                
        plotWidget.canvas.draw()

        
    def plot2dError(self, X, Y, T, Opt='XY'):        
        # Systematic and Random Errors ------
        plotWidget = self.plotWidgetXY
        xLabel = 'X(mm)'
        yLabel = 'Y(mm)'
        Xo = self.RefPoint[0]
        Yo = self.RefPoint[1]
        
        if Opt == 'XY':
            plotWidget = self.plotWidgetXY
            xLabel = 'X(mm)'
            yLabel = 'Y(mm)'
            Xo = self.RefPoint[0]
            Yo = self.RefPoint[1]
            
        if Opt == 'YZ':
            plotWidget = self.plotWidgetYZ
            xLabel = 'Y(mm)'
            yLabel = 'Z(mm)'
            Xo = self.RefPoint[1]
            Yo = self.RefPoint[2]            
            
        if Opt == 'XZ':
            plotWidget = self.plotWidgetXZ
            xLabel = 'X(mm)'
            yLabel = 'Z(mm)'
            Xo = self.RefPoint[0]
            Yo = self.RefPoint[2]            
                
        Xbeam, Ybeam = [], []
        if self.checkBoxBeamOn.isChecked() or self.checkBoxBeamEach.isChecked():
            Ts = NUMPY.asarray(self.TxTimeStartedList)
            Te = NUMPY.asarray(self.TxTimeStoppedList)    
            # print 'I, Ts Te', T, Ts, Te            
            i, j = 0, 0
            for i in range(len(T)):
                for j in range(len(Ts)):
                    if T[i] >= Ts[j] and T[i] <= Te[j]:
                        # print Opt,': ', Ts[j], T[i], Te[j], X[i], Y[i]
                        Xbeam.append(X[i])
                        Ybeam.append(Y[i])
                        break

            #print 'Xbeam Ybeam', Xbeam, Ybeam
            
        N = len(T)
                        
        plotWidget.canvas.ax.clear()   
        
        [Xmean, Xstd] = [Xo,0]
        [Ymean, Ystd] = [Yo,0]
        
        if self.checkBoxBeamOverall.isChecked():
            # print 'plot2dError:: ',Opt,' R: mean +/- std = ', NUMPY.mean(R), NUMPY.std(R)
            [Xmean, Xstd] = self.MeanStdOverall(T, X)
            [Ymean, Ystd] = self.MeanStdOverall(T, Y)
            # print 'plot2dError:: ',Opt+" [Xmean, Xstd] = ", [Xmean, Xstd], " [Ymean, Ystd] = ", [Ymean, Ystd]
            
            #print 'Time\tX\tY\t', self.T0
            #for i in range(len(X)):
            #    print T[i],'\t', X[i],'\t',Y[i]
        
        elif self.checkBoxBeamOn.isChecked() or self.checkBoxBeamEach.isChecked():
            xMeanStd = self.MeanStdEachBeam(T, X)  
            yMeanStd = self.MeanStdEachBeam(T, Y)  
            xMeanStdSumBeams = self.MeanStdSumBeams(xMeanStd)
            yMeanStdSumBeams = self.MeanStdSumBeams(yMeanStd)
    
            [Xmean, Xstd] = [xMeanStdSumBeams[0], xMeanStdSumBeams[1]]
            [Ymean, Ystd] = [yMeanStdSumBeams[0], yMeanStdSumBeams[1]]
            # print 'plot2dError:: ',Opt+" [Xmean+, Xstd+] = ", [Xmean, Xstd], " [YmeanA, YstdA] = ", [Ymean, Ystd]
        
        if self.checkBoxStdMean.isChecked():
            ellipse = Ellipse(xy=[Xmean, Ymean], width=Xstd*2, height=Ystd*2, angle = 0)
            ellipse.set_alpha(0.4)
            ellipse.set_facecolor('yellow')
            plotWidget.canvas.ax.add_patch(ellipse)
            # print 'plot2dError:: Xmean, Ymean, Xstd, Ystd = ', Xmean, Ymean, Xstd, Ystd
        
        
        
        [Xm, Xs] = [Xo,0]
        [Ym, Ys] = [Yo,0]
        if self.checkBoxBeamOverall.isChecked():
            # print 'plot2dError:: NUMPY.float64(Xo) = ', NUMPY.float64(Xo)
            # print 'plot2dError:: NUMPY.float64(Yo) = ', NUMPY.float64(Yo)
            [Xm, Xs] = self.MeanStdOverall(T, X, NUMPY.float64(Xo))
            [Ym, Ys] = self.MeanStdOverall(T, Y, NUMPY.float64(Yo))  
            # print 'plot2dError::Overall ',Opt+" [Xm, Xs] = ", [Xm, Xs], " [Ym, Ys] = ", [Ym, Ys]
              
        if self.checkBoxBeamOn.isChecked() or self.checkBoxBeamEach.isChecked():
            xMeanStd = self.MeanStdEachBeam(T, X, NUMPY.float64(Xo))  
            yMeanStd = self.MeanStdEachBeam(T, Y, NUMPY.float64(Yo))  
            xMeanStdSumBeams = self.MeanStdSumBeams(xMeanStd)
            yMeanStdSumBeams = self.MeanStdSumBeams(yMeanStd)
    
            [Xm, Xs] = [xMeanStdSumBeams[0], xMeanStdSumBeams[1]]
            [Ym, Ys] = [yMeanStdSumBeams[0], yMeanStdSumBeams[1]]
            # # print 'plot2dError::BeamOn ',Opt+" [Xm+, Xs+] = ", [Xm, Xs], " [Ym+, Ys+] = ", [Ym, Ys]
            
            X, Y = Xbeam, Ybeam
            
                
        if self.checkBoxStdRef.isChecked():
            ellipse0 = Ellipse(xy=[Xo, Yo], width=Xs*2, height=Ys*2, angle = 0)
            ellipse0.set_alpha(0.4)
            ellipse0.set_facecolor('green')
            plotWidget.canvas.ax.add_patch(ellipse0)
            # print 'plot2dError:: Xm, Ym, Xs, Ys = ', Xm, Ym, Xs, Ys

        
        A = []
        for i in range(N-1):
            A.append((T[i+1]-T[i])/(max(T)-min(T)))
                
            
        A = A/NUMPY.max(A)
        # print 'Alpha A = ', A
        # print 'N N(A)', N, len(A)
                
        if self.checkBoxAlpha.isChecked() and self.checkBoxBeamOverall.isChecked():
            for i in range(N-1):
                plotWidget.canvas.ax.scatter(X[i], Y[i], alpha=A[i])
        else:
            plotWidget.canvas.ax.plot(X, Y, 'o', color = 'blue',  ms=3, alpha=0.5) 
               
        #if self.checkBoxStdMean.isChecked():
        plotWidget.canvas.ax.plot([Xmean], [Ymean], 'o', color = 'red',   ms=5)
        
        #if self.checkBoxStdRef.isChecked():
        plotWidget.canvas.ax.plot([Xo],    [Yo],    'o', color = 'green', ms=5)
        
        if self.checkBoxAnnotate.isChecked():
            dist = NUMPY.sqrt((Xm - Xmean)*(Xm - Xmean)+(Ym - Ymean)*(Ym - Ymean))
            s = r'$\Delta$ = '+'{:.2f}'.format(dist)+' mm'
            plotWidget.canvas.ax.annotate(s, xy=(10, 25), xycoords='axes points', 
                    horizontalalignment='left', verticalalignment='bottom', fontsize=7)
            
            s = r'$\sigma_{Rx}$ = '+'{:.2f}'.format(Xstd)+', '+r'$\sigma_{Ry}$ = '+'{:.2f}'.format(Ystd)+' mm'
            plotWidget.canvas.ax.annotate(s, xy=(10, 15), xycoords='axes points', 
                    horizontalalignment='left', verticalalignment='bottom', fontsize=7)
            
            s = r'$\sigma_{Mx}$ = '+'{:.2f}'.format(Xs)+', '+r'$\sigma_{My}$ = '+'{:.2f}'.format(Ys)+' mm'
            plotWidget.canvas.ax.annotate(s, xy=(10, 5), xycoords='axes points', 
                    horizontalalignment='left', verticalalignment='bottom', fontsize=7)        
        

        if self.checkBoxWrite.isChecked():
            if Opt == 'XY': 
                self.Xmean, self.Ymean = Xmean, Ymean
                self.Xstd,  self.Ystd  = Xstd,  Ystd
                self.Xdist, self.Ydist = NUMPY.abs(Xmean-Xo), NUMPY.abs(Ymean-Yo)              
                #txtFileRst.write('\t'+Opt+'\t'+str(Xmean)+'\t'+str(Ymean)+'\t'\
                #                              +str(Xstd)+'\t'+str(Ystd)+'\t'\
                #                              +str(NUMPY.abs(Xmean-Xo))+'\t'+str(NUMPY.abs(Ymean-Yo)))
            if Opt == 'YZ': 
                self.Ymean, self.Zmean = Xmean, Ymean
                self.Ystd,  self.Zstd  = Xstd,  Ystd
                self.Ydist, self.Zdist = NUMPY.abs(Xmean-Xo), NUMPY.abs(Ymean-Yo) 
                #txtFileRst.write('\t'+Opt+'\t'+str(Xmean)+'\t'+str(Ymean)+'\t'\
                #                              +str(Xstd)+'\t'+str(Ystd)+'\t'\
                #                              +str(NUMPY.abs(Xmean-Xo))+'\t'+str(NUMPY.abs(Ymean-Yo)))
            if Opt == 'XZ':
                self.Xmean, self.Zmean = Xmean, Ymean
                self.Xstd,  self.Zstd  = Xstd,  Ystd
                self.Xdist, self.Zdist = NUMPY.abs(Xmean-Xo), NUMPY.abs(Ymean-Yo) 
                #txtFileRst.write('\t'+Opt+'\t'+str(Xmean)+'\t'+str(Ymean)+'\t'\
                #                              +str(Xstd)+'\t'+str(Ystd)+'\t'\
                #                              +str(NUMPY.abs(Xmean-Xo))+'\t'+str(NUMPY.abs(Ymean-Yo))+'\n')      
             
        
        labels_x = plotWidget.canvas.ax.get_xticklabels()
        labels_y = plotWidget.canvas.ax.get_yticklabels()
        for xlabel in labels_x:
            xlabel.set_fontsize(7)

        for ylabel in labels_y:
            ylabel.set_fontsize(7)
                              
        
        if self.checkBoxAspect.isChecked() and Opt=='XZ':
            # print 'XY lim = ', self.plotWidgetXY.canvas.ax.get_xlim(), self.plotWidgetXY.canvas.ax.get_ylim()
            # print 'YZ lim = ', self.plotWidgetYZ.canvas.ax.get_xlim(), self.plotWidgetYZ.canvas.ax.get_ylim()
            # print 'XZ lim = ', self.plotWidgetXZ.canvas.ax.get_xlim(), self.plotWidgetXZ.canvas.ax.get_ylim()
            xlimXY, ylimXY = self.plotWidgetXY.canvas.ax.get_xlim(), self.plotWidgetXY.canvas.ax.get_ylim()
            xlimYZ, ylimYZ = self.plotWidgetYZ.canvas.ax.get_xlim(), self.plotWidgetYZ.canvas.ax.get_ylim()
            xlimXZ, ylimXZ = self.plotWidgetXZ.canvas.ax.get_xlim(), self.plotWidgetXZ.canvas.ax.get_ylim()
            xlimRangeXY, ylimRangeXY = xlimXY[1] - xlimXY[0], ylimXY[1] - ylimXY[0]
            xlimRangeYZ, ylimRangeYZ = xlimYZ[1] - xlimYZ[0], ylimYZ[1] - ylimYZ[0]
            xlimRangeXZ, ylimRangeXZ = xlimXZ[1] - xlimXZ[0], ylimXZ[1] - ylimXZ[0]
            # print 'xlimRange XY YZ XZ = ', xlimRangeXY, xlimRangeYZ, xlimRangeXZ 
            # print 'ylimRange XY YZ XZ = ', ylimRangeXY, ylimRangeYZ, ylimRangeXZ
            xlimMaxRange = max(xlimRangeXY, xlimRangeYZ, xlimRangeXZ)
            ylimMaxRange = max(ylimRangeXY, ylimRangeYZ, ylimRangeXZ)
            maxLimSize = max(xlimMaxRange, ylimMaxRange)
            # print 'xlimMaxRange, ylimMaxRange, maxLimSize = ', xlimMaxRange, ylimMaxRange, maxLimSize
            
            # print 'maxLimSize - xlimRangeXY = ', maxLimSize - xlimRangeXY
            xlimXYnew = [xlimXY[0] - (maxLimSize - xlimRangeXY)*0.5, xlimXY[1] + (maxLimSize - xlimRangeXY)*0.5]
            ylimXYnew = [ylimXY[0] - (maxLimSize - ylimRangeXY)*0.5, ylimXY[1] + (maxLimSize - ylimRangeXY)*0.5]
            xlimYZnew = [xlimYZ[0] - (maxLimSize - xlimRangeYZ)*0.5, xlimYZ[1] + (maxLimSize - xlimRangeYZ)*0.5]
            ylimYZnew = [ylimYZ[0] - (maxLimSize - ylimRangeYZ)*0.5, ylimYZ[1] + (maxLimSize - ylimRangeYZ)*0.5]
            xlimXZnew = [xlimXZ[0] - (maxLimSize - xlimRangeXZ)*0.5, xlimXZ[1] + (maxLimSize - xlimRangeXZ)*0.5]
            ylimXZnew = [ylimXZ[0] - (maxLimSize - ylimRangeXZ)*0.5, ylimXZ[1] + (maxLimSize - ylimRangeXZ)*0.5]  
            
            self.plotWidgetXY.canvas.ax.set_xlim(xlimXYnew[0],xlimXYnew[1])
            self.plotWidgetXY.canvas.ax.set_ylim(ylimXYnew[0],ylimXYnew[1])
            self.plotWidgetYZ.canvas.ax.set_xlim(xlimYZnew[0],xlimYZnew[1])
            self.plotWidgetYZ.canvas.ax.set_ylim(ylimYZnew[0],ylimYZnew[1])
            self.plotWidgetXZ.canvas.ax.set_xlim(xlimXZnew[0],xlimXZnew[1])
            self.plotWidgetXZ.canvas.ax.set_ylim(ylimXZnew[0],ylimXZnew[1])

            #self.plotWidgetXY.canvas.ax.text(xlimXYnew[1]-maxLimSize/20.0,ylimXYnew[1]-maxLimSize/20.0, 'XY') 
            #self.plotWidgetYZ.canvas.ax.text(xlimYZnew[1]-maxLimSize/20.0,ylimYZnew[1]-maxLimSize/20.0, 'YZ')  
            #self.plotWidgetXZ.canvas.ax.text(xlimXZnew[1]-maxLimSize/20.0,ylimXZnew[1]-maxLimSize/20.0, 'XZ')   

                        
            self.plotWidgetXY.canvas.draw()
            self.plotWidgetYZ.canvas.draw()

        # print 'xLabel, yLabel = ', xLabel, yLabel
        plotWidget.canvas.ax.yaxis.set_label_coords(-0.12, 0.5) 
        # print     plotWidget.canvas.ax.get_position()
        plotWidget.canvas.ax.set_position([0.15,0.1,0.8,0.85])
        plotWidget.canvas.ax.set_xlabel(xLabel, fontsize = 7) 
        plotWidget.canvas.ax.set_ylabel(yLabel, fontsize = 7)

        plotWidget.canvas.draw()

        
    def plot3d(self, T, X, Y, Z): 
          
        N = len(T)
        A = []
        for i in range(N-1):
            # print i, T[i+1],T[i], T[i+1]-T[i], max(T)-min(T), (T[i+1]-T[i])/(max(T)-min(T))
            A.append((T[i+1]-T[i])/(max(T)-min(T)))
            
        A = A/NUMPY.max(A)
        # print A
                
        self.plot3dWidget.canvas.ax.clear()
        # self.plot3dWidget.canvas.ax.plot(X, Y, Z, marker='o', color = 'blue', ms=2)
        
        if self.checkBoxAlpha.isChecked():
            for i in range(N-1):
                self.plot3dWidget.canvas.ax.scatter(X[i], Y[i], Z[i], alpha=A[i])
        else:
            self.plot3dWidget.canvas.ax.scatter(X, Y, Z, alpha=0.5)
            
        Xo = self.RefPoint[0]
        Yo = self.RefPoint[1]
        Zo = self.RefPoint[2]
        print 'plot3d:: RefPoint = ', Xo, Yo, Zo
        self.plot3dWidget.canvas.ax.scatter([Xo], [Yo], [Zo], c='green', marker='*', edgecolor='green')
        

                        
        labels_x = self.plot3dWidget.canvas.ax.get_xticklabels()
        labels_y = self.plot3dWidget.canvas.ax.get_yticklabels()
        labels_z = self.plot3dWidget.canvas.ax.get_zticklabels()
        
        for xlabel in labels_x:
            xlabel.set_fontsize(7)

        for ylabel in labels_y:
            ylabel.set_fontsize(7)
            
        for zlabel in labels_z:
            zlabel.set_fontsize(7)
                        
            
        xLim = self.plot3dWidget.canvas.ax.get_xlim()
        yLim = self.plot3dWidget.canvas.ax.get_ylim()
        zLim = self.plot3dWidget.canvas.ax.get_zlim()

        self.plot3dWidget.canvas.ax.set_xlabel("x(mm)", fontsize = 9) 
        self.plot3dWidget.canvas.ax.set_ylabel("y(mm)", fontsize = 9)  
        self.plot3dWidget.canvas.ax.set_zlabel("z(mm)", fontsize = 9)  
        
        # self.plot3dWidget.canvas.ax.view_init(0,270)
        self.plot3dWidget.canvas.draw()
            
        self.plot3dWidget.canvas.ax.mouse_init()        
        

    def rotate3dXY(self):
        self.plot3dWidget.canvas.ax.view_init(90,-90)        
        self.plot3dWidget.canvas.draw()
        self.plot3dWidget.canvas.ax.mouse_init() 
        
    def rotate3dYZ(self):
        self.plot3dWidget.canvas.ax.view_init(0,0)        
        self.plot3dWidget.canvas.draw()
        self.plot3dWidget.canvas.ax.mouse_init() 
        
    def rotate3dXZ(self):
        self.plot3dWidget.canvas.ax.view_init(0,270)        
        self.plot3dWidget.canvas.draw()
        self.plot3dWidget.canvas.ax.mouse_init() 
    
    def reportShots(self):
        nRows = self.TableWidgetShot.rowCount()
        nColumns = self.TableWidgetShot.columnCount()
        
        self.TableLabelShot = [
            "ShotID", "Planned", "Elapsed", "Started", "Stopped", "Beam On",
            "Beam Off", "Duration", "R_Mean", "R_StdDev", "X_Mean", "X_StdDev",
            "Y_Mean", "Y_StdDev", "Z_Mean", "Z_StdDev"]
        
        iShotID   = self.TableLabelShot.index('ShotID')
        iPlanned  = self.TableLabelShot.index('Planned')
        iElapsed  = self.TableLabelShot.index('Elapsed')
        iDuration = self.TableLabelShot.index('Duration')
        iStarted  = self.TableLabelShot.index('Started')
        iStopped  = self.TableLabelShot.index('Stopped')
        iBeamOn   = self.TableLabelShot.index('Beam On')
        iBeamOff  = self.TableLabelShot.index('Beam Off')
        iMean     = self.TableLabelShot.index('R_Mean')
        iStdDev   = self.TableLabelShot.index('R_StdDev')
        iMeanX    = self.TableLabelShot.index('X_Mean')
        iStdDevX  = self.TableLabelShot.index('X_StdDev')
        iMeanY    = self.TableLabelShot.index('Y_Mean')
        iStdDevY  = self.TableLabelShot.index('Y_StdDev')
        iMeanZ    = self.TableLabelShot.index('Z_Mean')
        iStdDevZ  = self.TableLabelShot.index('Z_StdDev')
        
        
        HOME  = self.HomeLineEdit.text()
        gkDir = ''
        if self.GKDirListWidget.count() > 0:
            gkDir = self.GKDirListWidget.currentItem().text()

        gkLog = self.GKLogListWidget.currentItem().text()
        gkLogFile = gkLog
        if self.checkBoxWriteFile.isChecked():
            gkLogFile = str(DATETIME.date.today())
                        
        RptFile = str(HOME +'/'+ gkDir+'/'+ gkLogFile + '.rpt')
        txtFileRpt = open(RptFile, "w")
            
        # print 'ShotID\tTime\tX_Ref\tX_Dev\tY_Ref\tY_Dev\tZ_Ref\tZ_Dev'
        txtFileRpt.write('ShotID\tTime\tX_Ref\tX_Dev\tY_Ref\tY_Dev\tZ_Ref\tZ_Dev\n')
                    
        for i in range(0, nRows):
            if self.TableWidgetShot.item(i,iDuration) is not None:
                ShotID = self.TableWidgetShot.item(i,iShotID).text()
                BeamTime = float(self.TableWidgetShot.item(i,iDuration).text()) 
                RefX   = float(self.TableWidgetShot.item(i,iMeanX).text())
                DevX   = float(self.TableWidgetShot.item(i,iStdDevX).text())
                RefY   = float(self.TableWidgetShot.item(i,iMeanY).text())
                DevY   = float(self.TableWidgetShot.item(i,iStdDevY).text())
                RefZ   = float(self.TableWidgetShot.item(i,iMeanZ).text())
                DevZ   = float(self.TableWidgetShot.item(i,iStdDevZ).text())
        
                # print ShotID,'\t', BeamTime,'\t', RefX,'\t', DevX, RefY,'\t', DevY,'\t', RefZ,'\t', DevZ
                txtFileRpt.write(ShotID+'\t'+str(BeamTime)+'\t'+str(RefX)+'\t'+str(DevX)+'\t' \
                                 +str(RefY)+'\t'+str(DevY)+'\t'+str(RefZ)+'\t'+str(DevZ)+'\n')
                
        txtFileRpt.close()
                                
if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    
    V = app.desktop().screenGeometry()
    h = V.height()
    w = V.width()
    print("The screen resolution (width x height) is the following:")
    print(str(w) + " x " + str(h))
    
    plot = Plot_Widget()     
    plot.show()
    sys.exit(app.exec_())
