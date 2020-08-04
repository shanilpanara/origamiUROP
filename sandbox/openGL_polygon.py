import sys

from PyQt5 import QtWidgets, QtGui, QtCore, uic
from PyQt5.QtOpenGL import QGLWidget

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

vertices= (
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, -1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, -1, 1),
    (-1, 1, 1)
    )

edges = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,7),
    (6,3),
    (6,4),
    (6,7),
    (5,1),
    (5,4),
    (5,7)
    )

def cube(*args, **kwargs):
    def wrapper():
        glColor3f(1,0,0)
        glBegin(GL_LINES)
        for edge in edges:
            for vertex in edge:
                glVertex3fv(vertices[vertex])
        glEnd()
        gluPerspective(45, (args[0]/args[1]), 0.1, 50.0)
    return wrapper


class Window(QtWidgets.QMainWindow):
    def __init__(self, *args):
        super().__init__(*args)
        self.width = 640
        self.height = 480
        uic.loadUi('qt.ui', self)
        self.canvas = self.findChild(QGLWidget, 'canvas')
        self.setup_UI()

    def setup_UI(self):
        self.canvas.resize(self.width, self.height)
        self.canvas.initializeGL()
        self.canvas.resizeGL(self.width, self.height)
        self.canvas.paintGL = cube(self.width, self.height)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.canvas.update)
        timer.start(1000)
        return

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT)
        glRotatef(1, 3, 1, 1)
        glTranslatef(0.0, -0.0, -10.)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())