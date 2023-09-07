#include <GL/glut.h>
#include<cmath>
#include<vector>
#include <iostream>
#include "include/utils.h"

using namespace std;
int winWidth = 800;
int winHeight = 600;
static double boundary = 4;


static vector<Particle3D> particles;

void renderString(float x, float y, void* font, const char* string) {
    const char* c;
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, winWidth, 0.0, winHeight);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glRasterPos2f(x, y);
    for (c = string; *c != '\0'; c++) {
        glutBitmapCharacter(font, *c);
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}
void drawAxis()
{
    
    /*glColor3f(1.0, 0.0, 0.0);
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
    glPopMatrix();*/

    // 绘制坐标轴
// X轴: 红色
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(8.0, 0.0, 0.0);  
    glVertex3f(0.0, 4.0, 0.0);
    glVertex3f(4.0, 4.0, 0.0);
    glEnd();

    // Y轴: 绿色
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 8.0, 0.0);
    glVertex3f(4.0, 0.0, 0.0);
    glVertex3f(4.0, 4.0, 0.0);
    glEnd();

    // Z轴: 蓝色
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 8.0);
    glEnd();

    // 渲染2D标签
    renderString(winWidth * 0.1, winHeight *0.1, GLUT_BITMAP_9_BY_15, "X");
    renderString(winWidth *0.9, winHeight * 0.15, GLUT_BITMAP_9_BY_15, "Y");
    renderString(winWidth *0.5, winHeight *0.9, GLUT_BITMAP_9_BY_15, "Z");

    //坐标原点的球
    glPushMatrix(); // 保存当前的模型视图矩阵
    glColor3f(0.0, 0.0, 0.0);
    glutSolidSphere(0.2, 20, 20);
    glPopMatrix(); // 恢复之前的模型视图矩阵
}

//以米为单位，计算两个粒子之间的距离
static double distance(Particle3D& p1, Particle3D& p2)
{
    double res = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z,2);
    return sqrt(res) / 100;
}

static double pressure(Particle3D& p) {
    return 0.4 * (p.rho / RHO0 - 1);
}

//poly核函数用于密度计算
static double K_Poly6 = 315 / (64 * PI * pow(H, 9));
static double W_Poly6(double r)
{
    double res = K_Poly6 * pow((pow(H, 2) - pow(r, 2)), 3);
    return res;
}

//spiky核函数用于压力
static double K_Spiky = 45 / (PI * pow(H, 6));
static double W_Spiky(double r)
{
    double res = K_Spiky / r * pow(H - r, 2);
    return res;
}

//viscosity核用于粘度计算
static double K_Viscosity = 45 / (PI * pow(H, 6));
static double W_Viscosity(double r)
{
    double res = K_Viscosity * (H - r);
    return res;
}


static double getValueOfPressureAccerlate(Particle3D& i, Particle3D& j)
{
    double r = distance(i, j);
    double res = 0;
    if (r <= H)
    {
        //cout << "here,pressure is: " << (i.p + j.p) << endl;
        res = MASS * (i.p + j.p) / (2 * i.rho * j.rho) * W_Spiky(r);
    }
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
            particles[i].p = pressure(particles[i]);
        }
    }
}

static vector<double> getViscosityAcerlate(Particle3D& i, Particle3D& j)
{
    double r = distance(i, j);
    vector<double> res{ 0,0,0 };
    if (r <= H)
    {
        //cout << "W_Viscosity is: " << W_Viscosity(r) << endl;
        //cout << "multi rho is: " << i.rho * j.rho << endl;
        res[0] = MASS * MU * (j.vx - i.vx) / (i.rho * j.rho) * W_Viscosity(r);
        res[1] = MASS * MU * (j.vy - i.vy) / (i.rho * j.rho) * W_Viscosity(r);
        res[2] = MASS * MU * (j.vz - i.vz) / (i.rho * j.rho) * W_Viscosity(r);
    }
    return res;
}

static void calculateAcceleration(Particle3D& p)
{
    p.ax = 0.0;
    p.ay = 0.0;
    p.az = G;
    for (Particle3D& other : particles)
    {
        if (&p == &other) {
            continue;
        }
        double pressureAccerlate = getValueOfPressureAccerlate(p, other);
        //cout << "distance is: " << distance(p, other) << endl;
        //cout << "pressureAccerlate is: " << pressureAccerlate << endl;
        p.ax += pressureAccerlate * ((p.x - other.x)/100);
        p.ay += pressureAccerlate * ((p.y - other.y)/100);
        p.az += pressureAccerlate * ((p.z - other.z)/100);
        //cout << "pressure ax is: " << p.ax << " py is: " << p.ay <<" pz is: "<<p.az << endl;
        auto viscosityAccerlate = getViscosityAcerlate(p, other);
        p.ax += viscosityAccerlate[0];
        p.ay += viscosityAccerlate[1];
        p.az += viscosityAccerlate[2];
        //cout << "viscosity is: " << viscosityAccerlate[0] << "  " << viscosityAccerlate[1] <<"  "<<viscosityAccerlate[2] << endl;
        //cout << "final ax is: " << p.ax << " ay is: " << p.ay <<" az is: "<<p.az << endl;

    }
}

static void updateParticles()
{
    int i = 0;
    for (Particle3D& p : particles)
    {
        //cout << "here\n";
        //cout << "num is: " << i++ << endl;
        //cout << "position is: " << p.x << " " << p.y << " " << p.z << endl;
        //cout << "velocity is: " << p.vx << "  " << p.vy <<" "<<p.vz << endl;
        p.vx += p.ax * timeStep;
        p.vy += p.ay * timeStep;
        p.vz += p.az * timeStep;
        p.x += p.vx * timeStep * 100;
        p.y += p.vy * timeStep * 100;
        p.z += p.vz * timeStep * 100;
        if (p.x > boundary)
        {
            p.vx = -Vel_atten * p.vx;
            p.x = boundary;
        }
        if (p.x < 0)
        {
            p.vx = -Vel_atten * p.vx;
            p.x = 0;
        }
        if (p.y > boundary)
        {
            p.vy = -Vel_atten * p.vy;
            p.y = boundary;
        }
        if (p.y < 0)
        {
            p.vy = -Vel_atten * p.vy;
            p.y = 0;
        }
        if (p.z > boundary)
        {
            p.vz = -Vel_atten * p.vz;
            p.z = boundary;
        }
        if (p.z < 0)
        {
            p.vz = -Vel_atten * p.vz;
            p.z = 0;
        }
        //cout << p.ax << "  " << p.ay <<" "<<p.az << endl;
        //cout << "position is: " << p.x << " " << p.y <<" "<<p.z << endl;
        calculateAcceleration(p);
        //cout << endl;
    }
}

void initializeParticles3D(bool scale)
{
    /*particles.push_back({ 3, 3, 5, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0 });
    particles.push_back({ 3, 3, 4, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0 });*/
    double length = 0;
    if (scale)
        length = 0.4;
    for (double i = 2; i < 2.6+length; i+=0.2)
    {
        for (double j = 3; j < 3.6+ length; j+=0.2)
        {
            for (double k = 3; k < 3.6+length; k+=0.2)
            {
                Particle3D p;
                p.x = i;
                p.y = j;
                p.z = k;
                p.vx = -0.5;
                p.vy = -0.5;
                p.vz = 0;
                p.ax = 0.0;
                p.ay = 0.0;
                p.az = 0;
                p.rho = 0.0;
                p.p = 0.0;
                particles.push_back(p);
            }
        }
    }
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
    glutTimerFunc(1, timerCallback, 0);  // 每1ms重新设置定时器
}
void display3D() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawAxis();
    drawParticles();

    glutSwapBuffers();
}
int SPH_3D(int argc, char** argv,bool scale) {
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


    float dist = 12.0;
    float angle = 35.264;  // X和Z轴与观察平面的角度
    float isoX = dist * cos(angle * DEG_TO_RAD);
    float isoZ = dist * sin(angle * DEG_TO_RAD);
    float isoY = dist * sin(45 * DEG_TO_RAD);  // Y轴与观察平面的角度

    cout << "isoX is: " << isoX << " isoY is: " << isoY << " isoZ is: " << isoZ << endl;
    gluLookAt(isoX, isoY, isoZ, 0, 0, 2, 0, 0, 1);


    initializeParticles3D(scale);
    calculateDensityAndPressure();

    glutDisplayFunc(display3D);
    glutTimerFunc(1, timerCallback, 0);  // 设置初始定时器
    glutMainLoop();
    return 0;
}