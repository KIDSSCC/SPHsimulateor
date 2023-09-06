#include <GL/glut.h>
#include<cmath>
#include<vector>
#include <iostream>
#include "include/utils.h"

using namespace std;
int winWidth = 800;
int winHeight = 600;
double tmp = 5;

static vector<Particle3D> particles;


void drawAxis()
{
    glColor3f(1.0, 0.0, 0.0);
    glPushMatrix();
    glTranslatef(8, 0.0, 0.0);
    glutSolidCube(0.5);
    glPopMatrix();

    glColor3f(0.0, 1.0, 0.0);
    glPushMatrix();
    glTranslatef(0.0, 8, 0.0);
    glutSolidCube(0.5);
    glPopMatrix();

    glColor3f(0.0, 0.0, 1.0);
    glPushMatrix();
    glTranslatef(0.0, 0.0, 8);
    glutSolidCube(0.5);
    glPopMatrix();

    //坐标原点的球
    glPushMatrix(); // 保存当前的模型视图矩阵
    glColor3f(0.0, 0.0, 0.0);
    glutSolidSphere(0.5, 20, 20);
    glPopMatrix(); // 恢复之前的模型视图矩阵
}

//以米为单位，计算两个粒子之间的距离
static double distance(Particle3D& p1, Particle3D& p2)
{
    double res = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z,2);
    return sqrt(res) / 100;
}

//poly核函数用于密度计算
static double K_Poly6 = 315 / (64 * PI * pow(H, 9));
static double W_Poly6(double r)
{
    double res = K_Poly6 * pow((pow(H, 2) - pow(r, 2)), 3);
    return res;
}

static void calculateDensityAndPressure()
{
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i].rho = 0;
        for (int j = 0; j < particles.size(); j++)
        {
            double r = distance(particles[i], particles[j]);
            if (r < H)
            {
                particles[i].rho += MASS * W_Poly6(r);
            }
        }
    }
}

static void calculateAcceleration(Particle& p)
{

}

static void updateParticles()
{
    for (Particle3D& p : particles)
    {
        p.vx += p.ax * timeStep;
        p.vy += p.ay * timeStep;
        p.vz += p.az * timeStep;
        p.x += p.vx * timeStep * 100;
        p.y += p.vy * timeStep * 100;
        p.z += p.vz * timeStep * 100;
        if (p.x > 6)
        {
            p.vx = -Vel_atten * p.vx;
            p.x = 6;
        }
        if (p.x < 0)
        {
            p.vx = -Vel_atten * p.vx;
            p.x = 0;
        }
        if (p.y > 6)
        {
            p.vy = -Vel_atten * p.vy;
            p.y = 6;
        }
        if (p.y < 0)
        {
            p.vy = -Vel_atten * p.vy;
            p.y = 0;
        }
        if (p.z > 6)
        {
            p.vz = -Vel_atten * p.vz;
            p.z = 6;
        }
        if (p.z < 0)
        {
            p.vz = -Vel_atten * p.vz;
            p.z = 0;
        }
        calculateAcceleration(p);
    }
}

void initializeParticles3D()
{
    particles.push_back({ 3, 3, 5, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0 });
    particles.push_back({ 3, 3, 0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0 });
}

static void drawParticles()
{
    for (int i = 0; i < particles.size(); i++)
    {
        glPushMatrix(); // 保存当前的模型视图矩阵
        glColor3f(0.0, 0.0, 1.0);
        glTranslatef(particles[i].x, particles[i].y, particles[i].z);
        glutSolidSphere(0.07, 10, 10);
        glPopMatrix(); // 恢复之前的模型视图矩阵
    }
    
}

void timerCallback(int value) {
    updateParticles();
    calculateDensityAndPressure();
    timeElapsed += timeStep;

    glutPostRedisplay();
    glutTimerFunc(100, timerCallback, 0);  // 每1ms重新设置定时器
}
void display3D() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawAxis();
    drawParticles();

    glutSwapBuffers();
}
int SPH_3D(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(winWidth, winHeight);
    glutCreateWindow("3D Projection");
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.9, 0.9, 0.9, 1); // 设置背景颜色为浅灰色


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, (GLfloat)winWidth / winHeight, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    double camDistance = 20; // 调整这个值来放大或缩小视图
    double camX = camDistance * sin(45.0 * DEG_TO_RAD) * cos(35.264 * DEG_TO_RAD);
    double camY = camDistance * sin(45.0 * DEG_TO_RAD) * sin(35.264 * DEG_TO_RAD);
    double camZ = camDistance * cos(45.0 * DEG_TO_RAD);

    gluLookAt(camX, camY, camZ, 0, 0, 2, 0, 0, 1);


    initializeParticles3D();
    calculateDensityAndPressure();

    glutDisplayFunc(display3D);
    glutTimerFunc(100, timerCallback, 0);  // 设置初始定时器
    glutMainLoop();
    return 0;
}