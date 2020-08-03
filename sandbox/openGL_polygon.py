import sys

from PyQt5 import QtWidgets, QtGui, QtCore, uic

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

class Window(QtWidgets.QMainWindow):
    def __init__(self, *args):
        super().__init__(*args)
        uic.loadUi('qt.ui', self)
        self.canvas = self.findChild(QtWidgets.QOpenGLWidget, 'canvas')
        self.setup_UI()

    def setup_UI(self):
        self.canvas.initializeGL()
        self.canvas.resizeGL(1000,551)
        self.canvas.paintGL = self.paintGL
        #timer = QtCore.QTimer(self)
        #timer.timeout.connect(self.canvas.update)
        #timer.start(1000)
        return

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT)
        glColor3f(1,0,0)
        glBegin(GL_TRIANGLES)
        glVertex3f(-0.5,-0.5,0)
        glVertex3f(0.5,-0.5,0)
        glVertex3f(0.0,0.5,0)
        glEnd()

        gluPerspective(45, 1000/551, 0.1, 50.0)
        glTranslatef(0.0,0.0, -5)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())