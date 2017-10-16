import sys
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QWidget, QToolTip, QPushButton, QApplication, QVBoxLayout, QTabWidget, QLineEdit, QColorDialog, 
                             QSlider, QInputDialog, QFileDialog)
from PyQt5.QtGui import QFont, QColor   
from PyQt5.QtCore import QCoreApplication, Qt
from matplotlib import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import xml.etree.cElementTree as ET

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
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.mouseMoveEvent)
        self.canvas.mpl_connect('button_press_event', self.mousePrintEvent)
        #self.figure.mouseMoveEvent = self.mouseMoveEvent
        #self.figure.setMouseTracking(True)
        #self.canvas.viewport().installEventFilter(self)
        #self.setMouseTracking(True)

        self.ax = self.figure.add_subplot(111)  # create an axis
        self.ax.set_xlim([-100, 100])
        self.ax.set_ylim([-100, 100])

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

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        #self.setLayout(layout)
        self.canvas.draw()
        
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
 
    #def openFileNamesDialog(self):    
    #    options = QFileDialog.Options()
    #    options |= QFileDialog.DontUseNativeDialog
    #    files, _ = QFileDialog.getOpenFileNames(self,"QFileDialog.getOpenFileNames()", "","All Files (*);;Python Files (*.py)", options=options)
    #    if files:
    #        print(files)

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
            self.radius = int(elem.find("Radius").text)
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

        
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
