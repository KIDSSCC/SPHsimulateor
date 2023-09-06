#include <GL/glut.h>
#include <iostream>
#include "include/utils.h"

using namespace std;
void init() {
    glClearColor(1.0, 1.0, 1.0, 1.0);  // 设置背景色为黑色
    glEnable(GL_DEPTH_TEST);           // 启用深度测试，很重要，用于隐藏面消除
}
void display3D() {
    cout << "here\n";
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // 清除颜色和深度缓冲

    // 绘制球
    glPushMatrix(); // 保存当前的模型视图矩阵
    glColor3f(0.0, 0.0, 0.0);
    glutSolidSphere(1.0, 20, 20);
    glPopMatrix(); // 恢复之前的模型视图矩阵

    // 绘制三个立方体
    glPushMatrix();
    glColor3f(1.0, 0.0, 0.0);
    glTranslatef(2, 0, 0); // 设置立方体的位置
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

    glutSwapBuffers(); // 交换缓冲，用于双缓冲模式
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
    glutMainLoop(); // 开始GLUT主循环
    return 0;
}