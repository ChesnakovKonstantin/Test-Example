import sys
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QWidget, QToolTip, QPushButton, QApplication, QVBoxLayout, QTabWidget, QLineEdit, QColorDialog, 
                             QSlider, QInputDialog, QFileDialog)
from PyQt5.QtGui import QFont, QColor   
from PyQt5.QtCore import QCoreApplication, Qt
from matplotlib import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import xml.etree.cElementTree as ET
import threading

def dest(r1, r2, r3, t1, t2, t3):
    return ((r1- t1)**2+(r2- t2)**2+(r3- t3)**2)**(0.5)

def vectorfield(w, t, p):
        r11, v11, r12, v12, r13, v13, r21, v21, r22, v22, r23, v23, r31, v31, r32, v32, r33, v33 = w
        #r11 = r1[0]
        #r12 = r1[1]
        #r13 = r1[2]

        #v11 = v1[0]
        #v12 = v1[1]
        #v13 = v1[2]

        #r21 = r2[0]
        #r22 = r2[1]
        #r23 = r2[2]

        #v21 = v2[0]
        #v22 = v2[1]
        #v23 = v2[2]

        #r31 = r3[0]
        #r32 = r3[1]
        #r33 = r3[2]

        #v31 = v3[0]
        #v32 = v3[1]
        #v33 = v3[2]
        G, m1, m2, m3 = p

        # Create f = (x1',y1',x2',y2'):
        f = [v11, 
             G*m2*(r21 - r11)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r31-r11)/(dest(r31, r32, r33, r11, r12, r13)**3),
             v12, 
             G*m2*(r22 - r12)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r32-r12)/(dest(r31, r32, r33, r11, r12, r13)**3),
             v13, 
             G*m2*(r23 - r13)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r33-r13)/(dest(r31, r32, r33, r11, r12, r13)**3),
             v21, 
             G*m1*(r11 - r21)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r31-r21)/(dest(r31, r32, r33, r21, r22, r23)**3),
             v22, 
             G*m1*(r12 - r22)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r32-r22)/(dest(r31, r32, r33, r21, r22, r23)**3),
             v23, 
             G*m1*(r13 - r23)/(dest(r21, r22, r23, r11, r12, r13)**3)+ G*m3*(r33-r23)/(dest(r31, r32, r33, r21, r22, r23)**3),
             v31, 
             G*m1*(r11 - r31)/(dest(r31, r32, r33, r11, r12, r13)**3)+ G*m2*(r21-r31)/(dest(r31, r32, r33, r21, r22, r23)**3),
             v32, 
             G*m1*(r12 - r32)/(dest(r31, r32, r33, r11, r12, r13)**3)+ G*m2*(r22-r32)/(dest(r31, r32, r33, r21, r22, r23)**3),
             v33, 
             G*m1*(r13 - r33)/(dest(r31, r32, r33, r11, r12, r13)**3)+ G*m2*(r23-r33)/(dest(r31, r32, r33, r21, r22, r23)**3)]
        return f

class Circles():
    def __init__(self):
        self.elements = []

    def Add(self, el):
        self.elements.append(el)


class Example(QWidget):
    
    def __init__(self):
        super().__init__()
        
        
        self.initUI()
        
    def initUI(self):

        self.radius = 10
        self.color = 'r'

        self.circles = Circles()

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = QWidget()	
        self.tab2 = QWidget()
        self.tabs.resize(300,200) 
 
        # Add tabs
        self.tabs.addTab(self.tab1,"Edit")
        self.tabs.addTab(self.tab2,"Model")
        #self.setMouseTracking(True)
        self.figure = plt.figure()
        self.figureEarth = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvasEarth = FigureCanvas(self.figureEarth)
        self.canvas.mpl_connect('motion_notify_event', self.mouseMoveEvent)
        self.canvas.mpl_connect('button_press_event', self.mousePrintEvent)
        
        self.button_group = QtWidgets.QButtonGroup() # Number group
        self.r0 = QtWidgets.QRadioButton("verlet")
        self.button_group.addButton(self.r0)
        self.r1 = QtWidgets.QRadioButton("scipy")
        self.button_group.addButton(self.r1)
        self.r2 = QtWidgets.QRadioButton("threading")
        self.button_group.addButton(self.r2)
        self.button_group.buttonClicked.connect(self.RadioButtonClicked)

        self.ax = self.figure.add_subplot(111)  # create an axis
        self.ax.set_xlim([-100, 100])
        self.ax.set_ylim([-100, 100])

        self.axEarth = self.figureEarth.add_subplot(111)  # create an axis
        self.axEarth.set_xlim([-2.0* (10**11), 2.0* (10**11)])
        self.axEarth.set_ylim([-2.0* (10**11), 2.0* (10**11)])

        # Just some button connected to  method
        self.buttonPlus = QPushButton('+')
        self.buttonPlus.clicked.connect(self.IncreaseAxes)
        self.buttonPlus.resize(self.buttonPlus.sizeHint())
        self.buttonPlus.move(50, 50)

        self.buttonMinus = QPushButton('-')
        self.buttonMinus.clicked.connect(self.DecreaseAxes)
        self.buttonMinus.resize(self.buttonMinus.sizeHint())
        self.buttonMinus.move(70, 70)

        self.buttonSave = QPushButton('Save to xml file', self)
        self.buttonSave.move(10,10)
        self.buttonSave.clicked.connect(self.saveFileDialog)

        self.buttonOpen = QPushButton('Load from xml file', self)
        self.buttonOpen.clicked.connect(self.openFileDialog)
        self.buttonOpen.resize(self.buttonMinus.sizeHint())
        self.buttonOpen.move(70, 70)

        self.buttonColor = QPushButton('Open color dialog', self)
        self.buttonColor.setToolTip('Opens color dialog')
        self.buttonColor.move(10,10)
        self.buttonColor.clicked.connect(self.openColorDialog)

        self.textboxX = QLineEdit(self)
        self.textboxX.move(20, 20)
        self.textboxX.resize(120,40)
        #self.textboxX.setText('0')

        self.textboxY = QLineEdit(self)
        self.textboxY.move(220, 20)
        self.textboxY.resize(120,40)
        #self.textboxY.setText('0')

        self.sld = QSlider(Qt.Horizontal, self)
        self.sld.setFocusPolicy(Qt.NoFocus)
        self.sld.setGeometry(30, 40, 100, 30)
        self.sld.setRange(1.0, 100.0)
        self.sld.setValue(10.0)
        self.sld.valueChanged.connect(self.changeSliderRadius)

        self.textboxSld = QLineEdit(self)
        self.textboxSld.move(420, 20)
        self.textboxSld.resize(120,40)
        self.textboxSld.setText("10")
        self.textboxSld.textChanged.connect(self.changeTextRadius)

        # set the layout
        self.layout = QVBoxLayout(self)

        self.tab1.layout = QVBoxLayout(self)
        self.tab1.layout.addWidget(self.textboxX)
        self.tab1.layout.addWidget(self.textboxY)
        self.tab1.layout.addWidget(self.buttonColor)
        self.tab1.layout.addWidget(self.sld)
        self.tab1.layout.addWidget(self.textboxSld)
        self.tab1.layout.addWidget(self.canvas)
        self.tab1.layout.addWidget(self.buttonPlus)
        self.tab1.layout.addWidget(self.buttonMinus)
        self.tab1.layout.addWidget(self.buttonOpen)
        self.tab1.layout.addWidget(self.buttonSave)
        self.tab1.setLayout(self.tab1.layout)

        self.tab2.layout = QVBoxLayout(self)
        self.tab2.layout.addWidget(self.r0)
        self.tab2.layout.addWidget(self.r1)
        self.tab2.layout.addWidget(self.r2)
        self.tab2.layout.addWidget(self.canvasEarth)
        self.tab2.setLayout(self.tab2.layout)

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        #self.setLayout(layout)
        self.canvas.draw()
        self.canvasEarth.draw()
        
        #btn = QPushButton('Button', self)
        #btn.setToolTip('This is a <b>QuitButton</b> widget')
        #btn.clicked.connect(self.changeButtonName)
        #btn.resize(btn.sizeHint())
        #btn.move(50, 50)       
        
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('Circles')    
        self.show()

    def changeSliderRadius(self, value):
        self.radius = value
        print(value)
        self.textboxSld.setText(str(self.radius))

    def changeTextRadius(self, value):
        if value != '' and self.IsFloat(value):
            radius = float(value)
            self.radius = radius
            print(value)
            if radius < self.sld.minimum():
                self.textboxSld.setText(str(self.sld.minimum()))
                self.radius = self.sld.minimum()
            elif radius > self.sld.maximum():
                self.textboxSld.setText(str(self.sld.maximum()))
                self.radius = self.sld.maximum()
            self.sld.setValue(int(radius))

    def openColorDialog(self):
        color = QColorDialog.getColor()
        self.color = color.name()
        if color.isValid():
            print(color.name())

    def mouseMoveEvent(self, e):
    
        print('mouseEvent')
        #print(e)
        if (e.inaxes):
            x = e.xdata
            y = e.ydata
            self.textboxX.setText("{0}".format(x))
            self.textboxY.setText("{0}".format(y))

    def changeButtonName(self):
       
       self.setWindowTitle(Circle.CircleName());

    def IncreaseAxes(self):
        xlim = self.ax.get_xlim()
        self.ax.set_xlim(np.multiply(xlim, 1.5))
        ylim = self.ax.get_ylim()
        self.ax.set_ylim(np.multiply(ylim, 1.5))
        self.canvas.draw()
        self.sld.setMaximum(self.sld.maximum()*1.5)
        self.sld.setMinimum(self.sld.minimum()*1.5)

    def DecreaseAxes(self):
        xlim = self.ax.get_xlim()
        self.ax.set_xlim(np.divide(xlim, 1.5))
        ylim = self.ax.get_ylim()
        self.ax.set_ylim(np.divide(ylim, 1.5))
        self.canvas.draw()
        self.sld.setMaximum(self.sld.maximum()/1.5)
        self.sld.setMinimum(self.sld.minimum()/1.5)

    def mousePrintEvent(self, event):
        #if event.button() == QtCore.Qt.LeftButton:
        #    print("Press!")
        #super(GraphicsView, self).mousePressEvent(event)
        print('circle')
        #Circle.Draw(0, 0, self.radius, self.color, self.figure)
        circle1 = MyCircle(event.xdata, event.ydata, self.radius, self.color)
        self.ax.add_artist(circle1)
        self.canvas.draw()
        self.circles.Add(circle1)

    def CreateXMLRoot(self):
        root = ET.Element("root")
        doc = ET.SubElement(root, "doc")

        fig = ET.SubElement(doc, "figure")
        ET.SubElement(fig, "X").text = str(self.ax.get_xlim()[1])
        ET.SubElement(fig, "Y").text = str(self.ax.get_ylim()[1])
        ET.SubElement(fig, "Color").text = self.color

        slider = ET.SubElement(doc, "slider")
        ET.SubElement(slider, "Max").text = str(self.sld.maximum())
        ET.SubElement(slider, "Min").text = str(self.sld.minimum())
        ET.SubElement(slider, "Radius").text = str(self.radius)

        circles = ET.SubElement(doc, "circles")
        i = 0
        for element in self.circles.elements:
            el = ET.SubElement(circles, "circle{0}".format(i))
            ET.SubElement(el, "X").text = str(element.x)
            ET.SubElement(el, "Y").text = str(element.y)
            ET.SubElement(el, "radius").text = str(element.radius)
            ET.SubElement(el, "color").text = str(element.color)
            i += 1
        return root

    def saveFileDialog(self): 
        tree = ET.ElementTree(self.CreateXMLRoot())
        
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Save File as", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            print(fileName)
        tree.write(fileName)
 
    def ParseXMLFile(self, filename):
        tree = ET.ElementTree(file=filename)
        for elem in tree.iter(tag = "figure"):
            x = float(elem.find("X").text)
            y = float(elem.find("Y").text)
            self.color = elem.find("Color").text
            self.ax.set_xlim([-x, x])
            self.ax.set_ylim([-y, y])
        for elem in tree.iter(tag = "slider"):
            max = float(elem.find("Max").text)
            min = float(elem.find("Min").text)
            self.sld.setRange(min, max)
            self.radius = int(float(elem.find("Radius").text))
            self.sld.setValue(self.radius)
            self.textboxSld.setText(str(self.radius))
        self.ax.clear()
        for elem in tree.iter(tag = "circles"):
            for el in elem:
                x1 = float(el.find("X").text)
                y1 = float(el.find("Y").text)
                r = float(el.find("radius").text)
                col = el.find("color").text
                print('circle')
                circle1 = MyCircle(x1, y1, r, col)
                self.ax.add_artist(circle1)
                self.canvas.draw()
                self.circles.Add(circle1)
 
    def openFileDialog(self):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Open File","","All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            print(fileName)
            self.ParseXMLFile(fileName)

    def RadioButtonClicked(self, button):
        MEarth = 5.9724 * (10**24)
        MMoon = 7.34767309 * (10**22)
        MSun = 1988500 * (10**24)
        if (self.r0.isChecked()):
            N = 400000
            dt = 120
            Earth = Cosmic(N, MEarth, 0, -1.496*(10**11), 29.783*(10**3), 0, 0, 0, [])
            Moon = Cosmic(N, MMoon, 0, -1.496*(10**11) - 384.467*(10**6), 29.783*(10**3) + 1022, 0, 0, 0, [])
            Sun = Cosmic(N, MSun, 0, 0, 0, 0, 0, 0, [])
            cosmics = [Sun, Earth, Moon]
            for obj in cosmics:
                for interactionObj in cosmics:
                    if (not interactionObj is obj):
                        obj.Interactions.append((interactionObj.M, interactionObj.R))
            verlet = Verlet(N, dt, cosmics)
            verlet.VerletMain()
            self.axEarth.plot(Sun.R[0,:], Sun.R[1, :], color = 'yellow')
            self.axEarth.plot(Earth.R[0,:], Earth.R[1, :], color = 'blue')
            self.axEarth.plot(Moon.R[0,:], Moon.R[1, :], color = 'gray')
            self.canvasEarth.draw()
            print(self.button_group.checkedId())
            print(self.button_group.checkedButton().text())

        if (self.r1.isChecked()):
            p = [6.67408 * (10**(-11)), 1988500 * (10**24), 5.9724 * (10**24), 7.34767309 * (10**22)]
            w0 = [0, 0, 0, 0, 0, 0, 
                  0, 29.783*(10**3), -1.496*(10**11), 0, 0, 0,
                  0, 29.783*(10**3) + 1022, -1.496*(10**11) - 384.467*(10**6), 0, 0, 0]
            t = [120*float(i) for i in range(400000)]
            wsol = odeint(vectorfield, w0, t, args=(p,))
            xSun = wsol[:, 0]
            ySun = wsol[:, 2]
            xEarth = wsol[:, 6]
            yEarth = wsol[:, 8]
            xMoon = wsol[:, 12]
            yMoon = wsol[:, 14]
            print(5)
            self.axEarth.plot(xSun, ySun, color = 'yellow')
            self.axEarth.plot(xEarth, yEarth, color = 'blue')
            self.axEarth.plot(xMoon, yMoon, color = 'gray')
            self.canvasEarth.draw()

        if (self.r2.isChecked()):
            N = 4000
            dt = 120
            Earth = Cosmic(N, MEarth, 0, -1.496*(10**11), 29.783*(10**3), 0, 0, 0, [])
            Moon = Cosmic(N, MMoon, 0, -1.496*(10**11) - 384.467*(10**6), 29.783*(10**3) + 1022, 0, 0, 0, [])
            Sun = Cosmic(N, MSun, 0, 0, 0, 0, 0, 0, [])
            cosmics = [Sun, Earth, Moon]
            for obj in cosmics:
                for interactionObj in cosmics:
                    if (not interactionObj is obj):
                        obj.Interactions.append((interactionObj.M, interactionObj.R))
            verlet = VerletThreads(N, dt, cosmics)
            verlet.VerletMain()
            self.axEarth.plot(Sun.R[0,:], Sun.R[1, :], color = 'yellow')
            self.axEarth.plot(Earth.R[0,:], Earth.R[1, :], color = 'blue')
            self.axEarth.plot(Moon.R[0,:], Moon.R[1, :], color = 'gray')
            self.canvasEarth.draw()

            #verlet = VerletThreads(400000,120)
            #verlet.REarth[0,0] = 0
            #verlet.REarth[1,0] = -1.496*(10**11)
            #verlet.VEarth[0,0] = 29.783*(10**3)

            #verlet.RMoon[1,0] = -1.496*(10**11) - 384.467*(10**6)
            #verlet.VMoon[0,0] = 29.783*(10**3) + 1022
            #verlet.VerletMain()
            #print(5)
            #self.axEarth.plot(verlet.REarth[0,:], verlet.REarth[1, :], color = 'blue')
            #self.axEarth.plot(verlet.RMoon[0,:], verlet.RMoon[1, :], color = 'gray')
            #print(6)
            #self.canvasEarth.draw()

    def IsFloat(self, value):
        try:
            float(value)
            return True
        except:
            return False

class MyCircle(Circle):

    def __init__(self, x, y, radius, color):
        super().__init__((x, y), radius, color = color)
        self.x = x
        self.y = y
        self.radius = radius
        self.color = color

    def CircleName():
        return "circle"

    def Draw(x, y, radius, color, plt):
        circle1 = plt.Circle((x, y), radius, color = color)
        plt.ax
        ax = plt.gca()
        ax.add_artist(circle1)

class Cosmic:
    def __init__(self, n, mass, RxInit, RyInit, VxInit, VyInit, AxInit, AyInit, interactions):
        self.N = n
        self.M = mass
        self.R = np.zeros(shape = (3, self.N+1))
        self.V = np.zeros(shape = (3, self.N+1))
        self.A = np.zeros(shape = (3, self.N+1))
        self.R[0,0] = RxInit
        self.R[1,0] = RyInit
        self.V[0,0] = VxInit
        self.V[1,0] = VyInit
        self.A[0,0] = AxInit
        self.A[1,0] = AyInit
        self.Interactions = interactions

class Verlet:
    def __init__(self, Num, Dt, objects):
        self.N = Num
        self.dt = Dt
        self.G = 6.67408 * (10**(-11))
        self.Objects = objects

    def Acceleration(self, coordFirst, coordSecond, mass):
        if (np.linalg.norm(coordSecond - coordFirst) < 0.0000001):
            return 0
        return self.G*mass*(coordSecond - coordFirst)/(np.linalg.norm(coordSecond - coordFirst)**3)

    def VerletStep(self, R, V, A, i, kwargs):
        for (mass, r) in kwargs:
           A[:, i] += self.Acceleration(R[:, i], r[:, i], mass)
        R[:, i+1] = R[:, i] + V[:, i]*self.dt + (A[:, i]*(self.dt**2)*0.5)
        V[:, i+1] = V[:, i] + A[:, i]*self.dt

    def VerletMain(self):
        for i in range(self.N):
            for obj in self.Objects:
                interactions = []
                #for interactionObj in self.Objects:
                #    if (not interactionObj is obj):
                #        interactions.append((interactionObj.M, interactionObj.R))
                self.VerletStep(obj.R, obj.V, obj.A, i, obj.Interactions)
                print(i)

class VerletThreads:
    def __init__(self, Num, Dt, objects):
        self.N = Num
        self.dt = Dt
        self.G = 6.67408 * (10**(-11))
        self.Objects = objects

    def Acceleration(self, coordFirst, coordSecond, mass):
        if (np.linalg.norm(coordSecond - coordFirst) < 0.0000001):
            return 0
        return self.G*mass*(coordSecond - coordFirst)/(np.linalg.norm(coordSecond - coordFirst)**3)

    def VerletStep(self, R, V, A, eventMy, eventOther, kwargs):
        for i in range(self.N):
            for (mass, r) in kwargs:
                A[:, i] += self.Acceleration(R[:, i], r[:, i], mass)
            R[:, i+1] = R[:, i] + V[:, i]*self.dt + (A[:, i]*(self.dt**2)*0.5)
            V[:, i+1] = V[:, i] + A[:, i]*self.dt
            print(i)
            eventMy.set()
            eventOther.wait()
            eventOther.clear()
        eventMy.set()

    def MainThreadOperation(self, eventsMy, eventsOther):
        for i in range(self.N):
            #print("Main1")
            for evt in eventsMy:
                evt.wait()
                evt.clear()
            #print("Main2")
            for evt in eventsOther:
                evt.set()
            print("Main3")

    def VerletMain(self):
        eventMy = []
        eventTheir = []
        threads = []
        threadMain = threading.Thread(target = self.MainThreadOperation, args = (eventMy, eventTheir))

        for obj in self.Objects:
            interactions = []
            #for interactionObj in self.Objects:
            #    if (not interactionObj is obj):
            #       interactions.append((interactionObj.M, interactionObj.R))
            e = threading.Event()
            e.clear()
            eventTheir.append(e)
            ev = threading.Event()
            ev.clear()
            eventMy.append(ev)
            threads.append(threading.Thread(target = self.VerletStep, args = (obj.R, obj.V, obj.A, ev, e, obj.Interactions)))
            print(-1)
        for thrd in threads:
            thrd.start()
        threadMain.start()
        for thrd in threads:
            thrd.join()
        #threadMain.join()


        #event12 = threading.Event()
        #event12.clear()
        #event21 = threading.Event()
        #event21.clear()
        #t1 = threading.Thread(name = "first", target = self.VerletStep, args = (self.REarth, self.VEarth, self.AEarth, 1, event12, event21, {self.mSun : self.RSun, self.mMoon : self.RMoon}))
        #t2 = threading.Thread(name = "second", target = self.VerletStep, args = (self.RMoon, self.VMoon, self.AMoon, 2, event21, event12, {self.mSun : self.RSun, self.mEarth : self.REarth}))
        #t1.start()
        #t2.start()
        #t1.join()
        #t2.join()
        
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
