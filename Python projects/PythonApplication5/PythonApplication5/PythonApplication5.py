from OpenGL.GL import * 
from OpenGL.GLU import * 
from OpenGL.GLUT import * 
import sys
import math

import pyglet

import pywavefront

class renderParam( object ):
	def __init__( self ):
		self.initialColor = [1, 1, 1]
		self.drawColor = self.initialColor
		self.tVec = [0, 0, 0]
		self.mouseButton = None

	#def reset( self ):
	#	self.drawColor = self.initialColor
	#	setList( self.tVec, 0 )
	#	self.mouseButton = None

rP = renderParam()
#print(os.path.abspath(__file__))
#meshes = pywavefront.Wavefront('earth.obj')
#window = pyglet.window.Window(1024, 720, caption = 'Demo', resizable = True)

oldMousePos = [ 0, 0 ]
def mouseButton( button, mode, x, y ):
	global rP, oldMousePos
	if(mode == GLUT_DOWN):
		rP.mouseButton = button
	else:
		rP.mouseButton = None
	oldMousePos[0], oldMousePos[1] = x, y
	glutPostRedisplay()

def mouseMotion(x, y):
    global rP, oldMousePos
    deltaX = x - oldMousePos[0]
    deltaY = y - oldMousePos[1]
    print('button')
    print(rP.mouseButton)
    if(rP.mouseButton == GLUT_LEFT_BUTTON):
        factor = 0.001

        rP.tVec[0] = deltaX*factor
        rP.tVec[1] = -deltaY*factor
        print('left button')
        oldMousePos[0], oldMousePos[1] = x, y
        glTranslatef(rP.tVec[0], rP.tVec[1], rP.tVec[2])
        glutPostRedisplay()

step = 1

def specialkeys(key, x, y):
    # Сообщаем о необходимости использовать глобального массива pointcolor
    global pointcolor
    # Обработчики специальных клавиш
    if key == GLUT_KEY_UP:          # Клавиша вверх
        glRotatef(5, 1, 0, 0)  
        #glTranslatef(0.1, 0, 0)
   # Вращаем на 5 градусов по оси X
    if key == GLUT_KEY_DOWN:        # Клавиша вниз
        glRotatef(-5, 1, 0, 0)      # Вращаем на -5 градусов по оси X
    if key == GLUT_KEY_LEFT:        # Клавиша влево
        glRotatef(5, 0, 1, 0)       # Вращаем на 5 градусов по оси Y
    if key == GLUT_KEY_RIGHT:       # Клавиша вправо
        glRotatef(-5, 0, 1, 0)      # Вращаем на -5 градусов по оси Y
    if key == GLUT_KEY_END:         # Клавиша END
        # Заполняем массив pointcolor случайными числами в диапазоне 0-1
        pointcolor = [[random(), random(), random()], [random(), random(), random()], [random(), random(), random()]]


# Процедура подготовки шейдера (тип шейдера, текст шейдера)
def create_shader(shader_type, source):
    # Создаем пустой объект шейдера
    shader = glCreateShader(shader_type)
    # Привязываем текст шейдера к пустому объекту шейдера
    glShaderSource(shader, source)
    # Компилируем шейдер
    glCompileShader(shader)
    # Возвращаем созданный шейдер
    return shader


# Процедура перерисовки
def draw():
    glClear(GL_COLOR_BUFFER_BIT)                    # Очищаем экран и заливаем серым цветом
    glEnableClientState(GL_VERTEX_ARRAY)            # Включаем использование массива вершин
    glEnableClientState(GL_COLOR_ARRAY)             # Включаем использование массива цветов
    # Первый параметр - сколько используется координат на одну вершину
    # Второй параметр - определяем тип данных для каждой координаты вершины
    # Третий парметр - определяет смещение между вершинами в массиве
    # Если вершины идут одна за другой, то смещение 0
    # Четвертый параметр - указатель на первую координату первой вершины в массиве
    glVertexPointer(3, GL_FLOAT, 0, cubedata)
    # Указываем, где взять массив цветов:
    # Параметры аналогичны, но указывается массив цветов
    glColorPointer(3, GL_FLOAT, 0, cubecolor)
    # Рисуем данные массивов за один проход:
    # Первый параметр - какой тип примитивов использовать (треугольники, точки, линии и др.)
    # Второй параметр - начальный индекс в указанных массивах
    # Третий параметр - количество рисуемых объектов (в нашем случае это 3 вершины - 9 координат)
    glRotated(270, 1, 0, 0)
    glDrawArrays(GL_TRIANGLES, 0, 36 + 2160)
    glRotated(-270, 1, 0, 0)
    #Draw Cone
    #glBegin(GL_TRIANGLES)
    glDisableClientState(GL_VERTEX_ARRAY)           # Отключаем использование массива вершин
    glDisableClientState(GL_COLOR_ARRAY)            # Отключаем использование массива цветов
    glutSwapBuffers()                               # Выводим все нарисованное в памяти на экран
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)


# Здесь начинется выполнение программы
# Использовать двойную буферезацию и цвета в формате RGB (Красный Синий Зеленый)
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
glutMouseFunc(mouseButton)
glutMotionFunc(mouseMotion)
# Задаем серый цвет для очистки экрана
glClearColor(0.2, 0.2, 0.2, 1)
#Enable depth test
glEnable(GL_DEPTH_TEST)
#Accept fragment if it closer to the camera than the former one
glDepthFunc(GL_LESS)
# Создаем вершинный шейдер:
# Положение вершин не меняется
# Цвет вершины - такой же как и в массиве цветов
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
# Определяем массив вершин (три вершины по три координаты)
pointdata = [[0, 0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0], [-0.7, 1, -0.3]]

cubedata = [[-0.5, -0.5, -0.5],
    [-0.5,-0.5, 0.5],
    [-0.5, 0.5, 0.5],
    [0.5, 0.5,-0.5],
    [-0.5,-0.5,-0.5],
    [-0.5, 0.5,-0.5], 
    [0.5,-0.5, 0.5],
    [-0.5,-0.5,-0.5],
    [0.5,-0.5,-0.5],
    [0.5, 0.5,-0.5],
    [0.5,-0.5,-0.5],
    [-0.5,-0.5,-0.5],
    [-0.5,-0.5,-0.5],
    [-0.5, 0.5, 0.5],
    [-0.5, 0.5,-0.5],
    [0.5,-0.5, 0.5],
    [-0.5,-0.5, 0.5],
    [-0.5,-0.5,-0.5],
    [-0.5, 0.5, 0.5],
    [-0.5,-0.5, 0.5],
    [0.5,-0.5, 0.5],
    [0.5, 0.5, 0.5],
    [0.5,-0.5,-0.5],
    [0.5, 0.5,-0.5],
    [0.5,-0.5,-0.5],
    [0.5, 0.5, 0.5],
    [0.5,-0.5, 0.5],
    [0.5, 0.5, 0.5],
    [0.5, 0.5,-0.5],
    [-0.5, 0.5,-0.5],
    [0.5, 0.5, 0.5],
    [-0.5, 0.5,-0.5],
    [-0.5, 0.5, 0.5],
    [0.5, 0.5, 0.5],
    [-0.5, 0.5, 0.5],
    [0.5,-0.5, 0.5]]

cubecolor = [
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [1, 1, 0],
    [1, 1, 0],
    ]

for i in range(360):
    #glColor3f(1.0,1.0,0.0)
    #masconecolorlist.append([1.0,1.0,0.0])
    #masconevertslist.append([0, 0, 1])
    cubecolor.append([1.0,1.0,0.0])
    cubedata.append([0, 0, 1])
    #glVertex3f( 0, 0, 1.0)
    #glColor3f(1.0,1.0,0.0)
    #masconecolorlist.append([1.0,1.0,0.0])
    #masconevertslist.append([math.sin((i*math.pi)/180), math.cos((i*math.pi)/180), 0])
    cubecolor.append([1.0,1.0,0.0])
    cubedata.append([math.sin((i*math.pi)/180), math.cos((i*math.pi)/180), 0])
    #glVertex3f(math.sin((i*math.pi)/180), math.cos((i*math.pi)/180), 0)
    #glColor3f(1.0,1.0,0.0)
    #masconecolorlist.append([1.0,1.0,0.0])
    #masconevertslist.append([math.sin(((i+step)*math.pi)/180), math.cos(((i+step)*math.pi)/180), 0])
    cubecolor.append([1.0,1.0,0.0])
    cubedata.append([math.sin(((i+step)*math.pi)/180), math.cos(((i+step)*math.pi)/180), 0])
    #glVertex3f(math.sin(((i+step)*math.pi)/180), math.cos(((i+step)*math.pi)/180), 0)
#glEnd()

#glBegin(GL_TRIANGLES)
for i in range(360):
    #glColor3f(1.0,0.0,0.0)
    #glVertex3f( 0, 0, 0)
    cubecolor.append([1.0, 0.0, 0.0])
    cubedata.append([0, 0, 0])
    #glColor3f(1.0, 0.0, 0.0)
    #glVertex3f(math.cos((i*math.pi)/180), math.sin((i*math.pi)/180), 0)
    cubecolor.append([1.0,0.0,0.0])
    cubedata.append([math.cos((i*math.pi)/180), math.sin((i*math.pi)/180), 0])
    #glColor3f(1.0,0.0,0.0)
    #glVertex3f(math.cos(((i+step)*math.pi)/180), math.sin(((i+step)*math.pi)/180), 0)
    cubecolor.append([1.0,0.0,0.0])
    cubedata.append([math.cos(((i+step)*math.pi)/180), math.sin(((i+step)*math.pi)/180), 0])
#glEnd()

#cubedata = [[-0.5, -0.5, -0.5],
#    [-0.5,-0.5, 0.5],
#    [-0.5, 0.5, 0.5],
#    [0.5, 0.5,-0.5],
#    [-0.5,-0.5,-0.5],
#    [-0.5, 0.5,-0.5], 
#    [0.5,-0.5, 0.5],
#    [-0.5,-0.5,-0.5],
#    [0.5,-0.5,-0.5],
#    [0.5, 0.5,-0.5],
#    [0.5,-0.5,-0.5],
#    [-0.5,-0.5,-0.5],
#    [-0.5,-0.5,-0.5],
#    [-0.5, 0.5, 0.5],
#    [-0.5, 0.5,-0.5],
#    [0.5,-0.5, 0.5],
#    [-0.5,-0.5, 0.5],
#    [-0.5,-0.5,-0.5],
#    [-0.5, 0.5, 0.5],
#    [-0.5,-0.5, 0.5],
#    [0.5,-0.5, 0.5],
#    [0.5, 0.5, 0.5],
#    [0.5,-0.5,-0.5],
#    [0.5, 0.5,-0.5],
#    [0.5,-0.5,-0.5],
#    [0.5, 0.5, 0.5],
#    [0.5,-0.5, 0.5],
#    [0.5, 0.5, 0.5],
#    [0.5, 0.5,-0.5],
#    [-0.5, 0.5,-0.5],
#    [0.5, 0.5, 0.5],
#    [-0.5, 0.5,-0.5],
#    [-0.5, 0.5, 0.5],
#    [0.5, 0.5, 0.5],
#    [-0.5, 0.5, 0.5],
#    [0.5,-0.5, 0.5]]
#cubedata.append(masconevertslist)
#print(cubedata)
#print(masconevertslist)
## Определяем массив цветов (по одному цвету для каждой вершины)
#pointcolor = [[1, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 1]]
#cubecolor = [
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [0, 1, 0],
#    [0, 1, 0],
#    [0, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [0, 1, 0],
#    [0, 1, 0],
#    [0, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 0, 0],
#    [1, 0, 0],
#    [1, 0, 0],
#    [1, 0, 0],
#    [1, 0, 0],
#    [1, 0, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    [1, 1, 0],
#    ]
#cubecolor.append(masconecolorlist)
# Запускаем основой цикл
glutMainLoop()

