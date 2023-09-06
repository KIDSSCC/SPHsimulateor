#include <GL/glut.h>
#include <iostream>
#include "include/utils.h"

using namespace std;
void init() {
    glClearColor(1.0, 1.0, 1.0, 1.0);  // ���ñ���ɫΪ��ɫ
    glEnable(GL_DEPTH_TEST);           // ������Ȳ��ԣ�����Ҫ����������������
}
void display3D() {
    cout << "here\n";
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // �����ɫ����Ȼ���

    // ������
    glPushMatrix(); // ���浱ǰ��ģ����ͼ����
    glColor3f(0.0, 0.0, 0.0);
    glutSolidSphere(1.0, 20, 20);
    glPopMatrix(); // �ָ�֮ǰ��ģ����ͼ����

    // ��������������
    glPushMatrix();
    glColor3f(1.0, 0.0, 0.0);
    glTranslatef(2, 0, 0); // �����������λ��
    glutSolidCube(0.5);
    glPopMatrix();

    glPushMatrix();
    glColor3f(0.0, 1.0, 0.0);
    glTranslatef(0, 2, 0);
    glutSolidCube(0.5);
    glPopMatrix();

    glPushMatrix();
    glColor3f(0.0, 0.0, 1.0);
    glTranslatef(0, 0, 2);
    glutSolidCube(0.5);
    glPopMatrix();

    glutSwapBuffers(); // �������壬����˫����ģʽ
}

int SPH_3D(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("3D Scene");

    init();

    glutDisplayFunc(display3D);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 800.0 / 600.0, 1.0, 100.0);
    gluLookAt(5, 5, 5, 0, 0, 0, 0, 1, 0);
    glutMainLoop(); // ��ʼGLUT��ѭ��
    return 0;
}