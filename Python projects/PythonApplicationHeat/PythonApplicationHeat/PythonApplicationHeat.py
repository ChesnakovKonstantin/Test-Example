import math
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from random import random
import numpy as np
import random
import pywavefront

from OpenGL.raw.GLU import gluLookAt

import scipy.optimize as so
import scipy.integrate as si
import copy

def create_shader(shader_type, source):
    shader = glCreateShader(shader_type)
    glShaderSource(shader, source)
    glCompileShader(shader)
    return shader


def draw():
       
    glMatrixMode(GL_PROJECTION)
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_COLOR_ARRAY)
    glVertexPointer(3, GL_FLOAT, 0, pointdata)
    glColorPointer(3, GL_FLOAT, 0, pointcolor)
    glDrawArrays(GL_TRIANGLES, 0,len(pointdata)*3)#10262)#2196)
    glDisableClientState(GL_VERTEX_ARRAY)
    glDisableClientState(GL_COLOR_ARRAY)
    glutSwapBuffers()
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)


lastx=0
lasty=0

def MouseMotion (x, y):
    global lastx, lasty
    glTranslatef((x-lastx)/300,-(y-lasty)/300,0)
    lastx = x
    lasty = y
    glutPostRedisplay ()

def MouseRotate (x, y):
    global lastx, lasty, pointdata
    glRotatef((x-lastx)/3,0,1,0)
    glRotatef((y-lasty)/3,1,0,0)
    lastx = x
    lasty = y
    glutPostRedisplay ()
def specialkeys(key, x, y):
    global pointcolor
    if key == GLUT_KEY_UP:          
        glRotatef(5, 1, 0, 0)       
    if key == GLUT_KEY_DOWN:        
        glRotatef(-5, 1, 0, 0)      
    if key == GLUT_KEY_LEFT:        
        glRotatef(5, 0, 1, 0)       
    if key == GLUT_KEY_RIGHT:       
        glRotatef(-5, 0, 1, 0)
    if key == GLUT_KEY_END:
        pointcolor = [[random(), random(), random()], [random(), random(), random()], [random(), random(), random()]]
    if key == GLUT_KEY_F1:
        for j in range(0,2):
            for i in range(len(list_of_vel)):
                Col =Temp_to_color(temperature.T0[i])
                for j in range(diapazon[i], diapazon[i + 1]):
                    pointcolor[j] = [Col, Col, Col]
                #for j in range(diapazon[4], diapazon[5]):
                #    pointcolor[j] = [1, 1, 1]
            temperature.next_step()
        print(temperature.T0)


def parseFile(file_name):
    count = 0
    start = 0
    list_of_count_of_vel = []
    list_of_vel = []
    list_of_triang = []
    index = 0
    lst_vel = []
    lst_f = []
    total = 1
    count_of_v = 0

    for line in open(file_name, 'r'):
        values = line.split()
        if len(values) < 2:
            continue
        if(values[0]== '#' and values[1] == 'object' and count != 0):
            list_of_count_of_vel.append(count)
            list_of_vel.append(lst_vel)
            list_of_triang.append(lst_f)
            index = index + 1
            total = total + count_of_v
            count_of_v = 0
            count = 0
            lst_vel = []
            lst_f = []
        if (values[0] == '#' and values[1] == 'object' and count == 0):
            start = 1

        if(values[0] == 'f' and count == 0):
            start = 1

        if(start == 1 and values[0] == 'f'):
            count = count + 1
            lst_f.append([float(values[1])-total,float(values[2])-total,float(values[3])-total])

        if (start == 1 and values[0] == 'v'):
            lst_vel.append([float(values[1]), float(values[2]), float(values[3])])
            count_of_v = count_of_v + 1

    list_of_vel.append(lst_vel)
    list_of_triang.append(lst_f)
    list_of_count_of_vel.append(count)
    #print("list_of_count_of_vel=",list_of_count_of_vel)
    #print("list_of_triang=", list_of_triang)
    #print("list_of_vel=", list_of_vel)
    #print("list_of_count",sum(list_of_count_of_vel),3*sum(list_of_count_of_vel) )
    return  list_of_count_of_vel,list_of_triang,list_of_vel


def Form_Triangles(list_of_vel,list_of_tri):
    triangles = []
    diapazon = np.zeros(len(list_of_vel) + 1, dtype=np.int)
    for i in range(len(list_of_vel)):
        diapazon[i + 1] = diapazon[i] + len(list_of_tri[i])
        for el in list_of_tri[i]:
            triangle = np.array([list_of_vel[i][int(el[0])], list_of_vel[i][int(el[1])], list_of_vel[i][int(el[2])]])
            triangles.append(triangle)
    return np.array(triangles), diapazon


def Temp_to_color(temp):
    B=5
    D=5
    if temp<50:
        return [np.exp(-(temp-50)**2/B), np.exp(-(temp)**2/B),np.exp(-(temp+50)**2/B)]
    return [np.exp(-(temp-62)**2/B), 0, np.exp(-(temp-60)**2/D)]

def SquareTr(triangle):
    l1 = np.linalg.norm(triangle[0] - triangle[1])
    l2 = np.linalg.norm(triangle[0] - triangle[2])
    l3 = np.linalg.norm(triangle[2] - triangle[1])
    p=(l1+l2+l3)/2
    return (p*(p-l1)*(p-l2)*(p-l3))**0.5

def SquaresOfParts(triangles, diapazon):
    squares = []
    for i in range(len(diapazon) - 1):
        square = 0
        for j in range(diapazon[i], diapazon[i+1]):
            square += SquareTr(triangles[j])
        squares.append(square)
    return squares

if __name__ == '__main__':
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
    glutInitWindowSize(300, 300)
    glutInitWindowPosition(50, 50)
    glutInit(sys.argv)
    glutCreateWindow(b"Shaders!")
    glutDisplayFunc(draw)
    glutIdleFunc(draw)
    glutSpecialFunc(specialkeys)
    glutPassiveMotionFunc(MouseMotion)
    glutMotionFunc(MouseRotate)
    glClearColor(0.2, 0.2, 0.2, 1)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_MULTISAMPLE)
    vertex = create_shader(GL_VERTEX_SHADER, """
        varying vec4 vertex_color;

        void main() {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
            vertex_color = gl_Color;
        }
        """)
    fragment = create_shader(GL_FRAGMENT_SHADER, """
        varying vec4 vertex_color;

        void main() {
            gl_FragColor = vertex_color;
        }
        """)

    program = glCreateProgram()
    glAttachShader(program, vertex)
    glAttachShader(program, fragment)
    glLinkProgram(program)
    glUseProgram(program)
    
    listOfDotsNumber, list_of_tri, list_of_vel = parseFile('model2.obj')
    NumberOfParts = len(listOfDotsNumber)
    print("M=",NumberOfParts)
    triangles, diapazon = Form_Triangles(list_of_vel, list_of_tri)
    print("triangles", triangles)
    print("range", diapazon)
    squareOfParts = SquaresOfParts(triangles, diapazon)
    print("squares", squareOfParts)
    triangles /= (2 * triangles.max())
    pointcolorTri = np.zeros((diapazon[len(list_of_vel)], 3, 3))
    from random import random
    for i in range(0, len(listOfDotsNumber)):
        m = random()
        k = random()
        for j in range(diapazon[i], diapazon[i + 1]):
            pointcolorTri[j] = [k, 0.0, m]
        #pointcolorTri[0] = [1, 1, 1]
    #print("triangles=",triangles, "len=", len(triangles), " x ",len(triangles[0]),"=",len(triangles)*len(triangles[0]))
    #print("diapazon=", diapazon)

    pointdata=triangles
    pointcolor=pointcolorTri

    def func_solve(T, t, lambada, Q_R, c, eps, S):
        NumberOfParts = len(c)
        right_part = np.zeros(NumberOfParts)
        StephBolC = 5.67
        for i in range(NumberOfParts):
            for j in range(NumberOfParts):
                if i != j:
                    right_part[i] -= lambada[i, j] * S[i, j] * (T[i] - T[j])
            right_part[i] -= eps[i] * S[i, i] * StephBolC * (T[i] / 100) ** 4
            right_part[i] += Q_R[i](t)
            right_part[i] /= c[i]
        return right_part
    class TempComputer:
        def __init__(self, lambada, Q_R, c, eps, S, tau, NumberOfParts):
            self.lambada = lambada
            self.Q_R = Q_R
            self.c = c
            self.eps = eps
            self.S = S
            self.counter = 0
            print("NumberOfParts=",NumberOfParts,"c=",c)
            T = self.T0 = so.fsolve(func_solve, np.zeros(NumberOfParts), args=(0, lambada, Q_R, c, eps, S,))
            self.T = copy.copy(self.T0)
            self.tau = tau
            self.NumberOfParts = NumberOfParts

        def next_step(self):
            Tm = np.linspace((self.counter - 1) * self.tau, self.counter * self.tau, 2)
            print("Tm=",Tm)
            self.counter += 1
            self.T = si.odeint(func_solve, self.T0, Tm, args=(self.lambada, self.Q_R, self.c, self.eps, self.S,))
            self.T0 = copy.copy(self.T[1])
            return self.T[1]
    squareOfElem=np.array([[  36.1672021, 4*np.pi    ,  0        , 0          , 0],
                           [  4*np.pi   , 219.8340987, 4*np.pi   , 0          , 0],
                           [  0.        , 4*np.pi    , 12.3660143, np.pi      , 0],
                           [  0.        , 0.         , np.pi     , 99.4076902, np.pi],
                           [  0.        , 0.         ,0.         , np.pi      , 268]])

    squareOfElem=np.array([[  36.1672021, 4*np.pi   ,  0        , 0      , 0],
                           [  4*np.pi   , 99.4076907,  0        , 4*np.pi, 0],
                           [  0.        , 0         , 12.3660143, np.pi  , np.pi],
                           [  0.        , 4*np.pi   , np.pi     , 268    , 0],
                           [  0.        , 0.        , np.pi     , 0      , 219.8341097]])

    print("square=", squareOfElem)

    #for i in range(NumberOfParts):
    #    for j in range(i + 1, NumberOfParts):
    #        temp = (squareOfElem[i, j] + squareOfElem[j, i]) / 2
    #        squareOfElem[i, j] = temp
    #        squareOfElem[j, i] = temp
    #print ("squareOfElem=",squareOfElem)
    eps = [0.05, 0.05, 0.01, 0.1, 0.1]
    c = [520, 520, 840, 900, 900]

    lambada = np.zeros((NumberOfParts, NumberOfParts))
    lambada[0, 1] = 20
    lambada[1, 0] = 20
    lambada[1, 3] = 130
    lambada[3, 1] = 130
    lambada[2, 3] = 10.5
    lambada[3, 2] = 10.5
    lambada[2, 4] = 119
    lambada[4, 2] = 119

    Q_R = []
    for i in range(NumberOfParts):
        f = lambda t: [0]
        Q_R.append(f)
    A = 2
    Q_R[0] = lambda t: [A * (20 + 3 * np.cos(t / 4))]
    tau = 10 ** 2
    print("pointcolor[6]=",pointcolor[6])
    temperature = TempComputer(lambada, Q_R, c, eps, squareOfElem, tau, NumberOfParts)
    for i in range(len(list_of_vel)):
        Col = Temp_to_color(temperature.T0[i])
        print('len(list_of_vel) = ', len(list_of_vel))
        print("temperature.T0[i]=",temperature.T0[i])
        print("Col=",Col)

        for j in range(diapazon[i], diapazon[i + 1]):
            pointcolor[j] = [Col, Col, Col]
            
    print(3*12+361*2*3)
    print("len(pointdata)_all=",len(pointdata)*len(pointdata[0]))
    print("len(pointdata)=", len(pointdata))
    print("len(pointcolor)=", len(pointcolor))

    #global eye
    #eye = np.zeros(3)
    #global lat
    #lat = 0
    #global lon
    #lon = np.arctan2(0, -1)
    glutMainLoop()
