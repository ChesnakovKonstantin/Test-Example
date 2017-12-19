#!/usr/bin/env python
"""This script shows another example of using the PyWavefront module."""
# This example was created by intrepid94
import sys
sys.path.append('..')
import ctypes

import pyglet
from pyglet.gl import *

import pywavefront
import numpy as np
import random

from pywavefront import Wavefront
from OpenGL.GL import * 
from OpenGL.GLU import * 
from OpenGL.GLUT import * 

global N

lastx=0
lasty=0

def MouseMotion (x, y):
    global lastx, lasty
    glTranslatef(-(x-lastx)/300,(y-lasty)/300,0)
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


def create_shader(shader_type, source):
    # Создаем пустой объект шейдера
    shader = glCreateShader(shader_type)
    # Привязываем текст шейдера к пустому объекту шейдера
    glShaderSource(shader, source)
    # Компилируем шейдер
    glCompileShader(shader)
    # Возвращаем созданный шейдер
    return shader

def draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                    
    glEnableClientState(GL_VERTEX_ARRAY)            
    glEnableClientState(GL_COLOR_ARRAY)            
    glVertexPointer(3, GL_FLOAT, 0, pointdata)
    glColorPointer(3, GL_FLOAT, 0, pointcolor)
    glEnable(GL_DEPTH_TEST)
    glDrawArrays(GL_TRIANGLES,0,3*N)
    glDisableClientState(GL_VERTEX_ARRAY)           
    glDisableClientState(GL_COLOR_ARRAY)            
    glutSwapBuffers()     

glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
# Указываем начальный размер окна (ширина, высота)
glutInitWindowSize(300, 300)
# Указываем начальное
# положение окна относительно левого верхнего угла экрана
glutInitWindowPosition(50, 50)
# Инициализация OpenGl
glutInit(sys.argv)
# Создаем окно с заголовком "Shaders!"
glutCreateWindow(b"Shaders!")
# Определяем процедуру, отвечающую за перерисовку
glutDisplayFunc(draw)
# Определяем процедуру, выполняющуюся при "простое" программы
glutIdleFunc(draw)
# Определяем процедуру, отвечающую за обработку клавиш
glutSpecialFunc(specialkeys)
glutPassiveMotionFunc(MouseMotion)
glutMotionFunc(MouseRotate)
# Задаем серый цвет для очистки экрана
glClearColor(0.2, 0.2, 0.2, 1)
#Enable depth test
glEnable(GL_DEPTH_TEST)
#Accept fragment if it closer to the camera than the former one
glDepthFunc(GL_LESS)

vertex = create_shader(GL_VERTEX_SHADER, """
varying vec4 vertex_color;
            void main(){
                gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
                vertex_color = gl_Color;
            }""")
# Создаем фрагментный шейдер:
# Определяет цвет каждого фрагмента как "смешанный" цвет его вершин
fragment = create_shader(GL_FRAGMENT_SHADER, """
varying vec4 vertex_color;
            void main() {
                gl_FragColor = vertex_color;
}""")
# Создаем пустой объект шейдерной программы
program = glCreateProgram()
# Приcоединяем вершинный шейдер к программе
glAttachShader(program, vertex)
# Присоединяем фрагментный шейдер к программе
glAttachShader(program, fragment)
# "Собираем" шейдерную программу
glLinkProgram(program)
# Сообщаем OpenGL о необходимости использовать данную шейдерну программу при отрисовке объектов
glUseProgram(program)

name='earth.obj'
meshes = pywavefront.Wavefront(name)
ps=pywavefront.ObjParser(meshes,name)
ps.read_file(name)
pointdata2=ps.material.vertices
N=len(pointdata2)//24
pointdata=np.zeros((N,3,3))
pointcolor=np.zeros((N,3,3))
for i in range(0,N):
    for j in range(0,3):
        pointdata[i,j,0:3]=pointdata2[24*i+8*j+5:24*i+8*j+8]
pointdata/=2*pointdata.max()
for i in range(0,N//2):
    mas=[random.random(), random.random(), random.random()]
    pointcolor[2*i]=[mas, mas, mas]
    pointcolor[2*i+1]=pointcolor[2*i]
glutMainLoop()
