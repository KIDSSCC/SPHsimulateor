#include <GL/glut.h>
#include <cmath>
#include <vector>
#include<iostream>
#include "include/SPH2D.h"
#include "include/utils.h"
using namespace std;
// 常量的声明

std::vector<Particle*> particles;
double maxX = 1.0, maxY = 1.0;
double minX = 0.0, minY = 0.0;
double margin = 0;
bool gridOptimize = false;

//以米为单位，计算两个粒子之间的距离
double distance(Particle *p1, Particle *p2) {
    auto res = sqrt(pow((p1->x - p2->x) * winSize / 100, 2) + pow((p1->y - p2->y) * winSize / 100, 2));
    return res/100;
}


//poly核函数用于密度计算
double K_Poly6 = 4 / (PI * pow(H, 8));
double W_Poly6(double r)
{
    double res = K_Poly6 * pow((pow(H, 2) - pow(r, 2)), 3);
    return res;
}

//spiky核函数用于压力
double K_Spiky = 30 / (7 * PI * pow(H, 5));
double W_Spiky(double r)
{
    double res = K_Spiky /r*pow(H - r, 2);
    return res;
}

//viscosity核用于粘度计算
double K_Viscosity = 20 / (PI * pow(H, 5));
double W_Viscosity(double r)
{
    double res = K_Viscosity * (H - r);
    return res;
}

double getValueOfPressureAccerlate(Particle* i, Particle* j)
{
    double r = distance(i, j);
    //cout << "in getValueOfPressureAccerlate,distance is:" <<r<< endl;
    //cout << "相对位置：" << (i.y - j.y) * winSize / 100 << endl;
    double res = 0;
    if (r <= H)
    {
        //cout << "here,pressure is: " << (i.p + j.p) << endl;
        res = MASS * (i->p + j->p) / (2 * i->rho * j->rho) * W_Spiky(r);
    }
    return res;
}
vector<double> getViscosityAcerlate(Particle* i, Particle* j)
{
    double r = distance(i, j);
    vector<double> res{0,0};
    if (r <= H)
    {
        //cout << "W_Viscosity is: " << W_Viscosity(r) << endl;
        //cout << "multi rho is: " << i.rho * j.rho << endl;
        res[0] = MASS * MU * (j->vx - i->vx) / (i->rho * j->rho) * W_Viscosity(r);
        res[1] = MASS * MU * (j->vy - i->vy) / (i->rho * j->rho) * W_Viscosity(r);
    }
    return res;
}

//单个粒子处的压力
double pressure(Particle* p) {
    return 0.4 * (p->rho/RHO0 - 1);
}


// 计算加速度
void calculateAcceleration(Particle* p) {
    p->ax = 0.0;
    p->ay = G;

    if (gridOptimize)
    {
        //网格优化
        int curr_i = p->i;
        int curr_j = p->j;
        for (int i = max(curr_i - 1, 0); i <= min(curr_i + 1, 39); i++)
        {
            for (int j = max(curr_j - 1, 0); j <= min(curr_j + 1, 39); j++)
            {
                //从grid[i][j]中获取所有的粒子
                auto itea = grid[i][j]->next;
                for (int k = 0; k < grid[i][j]->num; k++)
                {
                    if (itea == p)
                        continue;

                    double pressureAccerlate = getValueOfPressureAccerlate(p, itea);
                    p->ax += pressureAccerlate * ((p->x - itea->x) * winSize / 100);
                    p->ay += pressureAccerlate * ((p->y - itea->y) * winSize / 100);
                    auto viscosityAccerlate = getViscosityAcerlate(p, itea);
                    p->ax += viscosityAccerlate[0];
                    p->ay += viscosityAccerlate[1];
                    itea = itea->next;
                }
            }
        }
    }
    else
    {
        for (Particle *other : particles) {
            if (p == other) {
                continue;
            }
            double pressureAccerlate = getValueOfPressureAccerlate(p, other);
            //cout << "distance is: " << distance(p, other) << endl;
            //cout << "pressureAccerlate is: " << pressureAccerlate << endl;
            p->ax += pressureAccerlate * ((p->x - other->x) * winSize / 100);
            p->ay += pressureAccerlate * ((p->y - other->y) * winSize / 100);
            //cout << "pressure ax is: " << p.ax << " py is: " << p.ay << endl;
            auto viscosityAccerlate = getViscosityAcerlate(p, other);
            p->ax += viscosityAccerlate[0];
            p->ay += viscosityAccerlate[1];
            //cout << "viscosity is: " << viscosityAccerlate[0] << "  " << viscosityAccerlate[1] << endl;
            //cout << "final ax is: " << p.ax << " ay is: " << p.ay << endl;
        }
    }
}

// 速度位置更新
void updateParticles() {
    int i = 0;
    for (Particle* p : particles) {
        
        //cout << "num is: " << i++ << endl;
        //cout << "position is: " << p.x << " " << p.y << endl;
        //cout << "velocity is: " << p.vx << "  " << p.vy << endl;
        p->vx += p->ax * timeStep;
        p->vy += p->ay * timeStep;
        p->x += p->vx * timeStep*10000/winSize;
        p->y += p->vy * timeStep*10000/winSize;
        if (p->x > maxX) {
            p->vx = -Vel_atten * p->vx;
            p->x = maxX;
        }
        if (p->x < minX) {
            p->vx = -Vel_atten * p->vx;
            p->x = minX;
        }
        if (p->y > maxY) {
            p->vy = -Vel_atten * p->vy;
            p->y = maxY;
        }
        if (p->y < minY) {
            p->vy = -Vel_atten * p->vy;
            p->y = minY;
        }
        //cout << p.ax << "  " << p.ay << endl;
        //cout << "position is: " << p.x << " " << p.y << endl;
        calculateAcceleration(p);
        //cout << endl;

        //更新网格情况
        joinInGrid(p);
    }
    //printGrid();
}


// 计算每个粒子处的密度和每个粒子自己的压力
void calculateDensityAndPressure() {
    for (Particle* par : particles)
    { 
        if (gridOptimize)
        {
            //将旧版的遍历全部更换为遍历局部网格内
            int curr_i = par->i;
            int curr_j = par->j;
            par->rho = 0;
            for (int i = max(curr_i - 1, 0); i <= min(curr_i + 1, 39); i++)
            {
                for (int j = max(curr_j - 1, 0); j <= min(curr_j + 1, 39); j++)
                {
                    //从grid[i][j]中获取所有的粒子
                    auto itea = grid[i][j]->next;
                    for (int k = 0; k < grid[i][j]->num; k++)
                    {
                        double r = distance(par, itea);
                        if (r < H)
                        {
                            par->rho += MASS * W_Poly6(r);
                        }
                        itea = itea->next;
                    }
                }
            }
            par->p = pressure(par);
        }
        else
        {
            par->rho = 0;
            for (Particle* other : particles)
            {
                double r = distance(par, other);
                //cout << "distance is: " << r << endl;
                if (r < H)
                {
                    //cout << "value is: " << W_Poly6(r) << endl;
                    par->rho += MASS * W_Poly6(r);
                }
            }
            par->p = pressure(par);
            //cout << "pressure is: " << par.p << endl;
            //cout << "rho is: " << par.rho << endl;
        }
    }

}

// Function to draw particles
void drawParticles() {
    int numSegment = 200;
    //绘制出的粒子半径为0.06cm，此时对应的光滑核函数的平滑半径为0.15cm，二者之间保持三倍关系
    float radius = 0.01;

    glBegin(GL_POINTS);
    for (Particle *p : particles) {
        glColor3f(0.0, 0.0, 1.0);
        for (int i = 0; i < numSegment; i++)
        {
            float angle = i * 2.0 * PI / numSegment;
            float x = p->x + radius * cos(angle);
            float y = p->y + radius * sin(angle);
            glVertex2f(x, y);
        }
    }
    glEnd();
}

// Function to display particles
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    drawParticles();
    glutSwapBuffers();
}

// Function to update particles and time elapsed
void update(int value) {
    updateParticles();
    calculateDensityAndPressure();
    timeElapsed += timeStep;
    glutPostRedisplay();
    glutTimerFunc(1, update, 0);
}

// Function to initialize OpenGL
void initOpenGL() {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(minX-margin, maxX+ margin, minY- margin, maxY+ margin);
}


// 粒子初始化
void initializeParticles(bool scale) {
    initGrid(40);
    /*particles.push_back(new Particle( 0.4, 0.2, 0, 0.0, 0.0, 0.0, 0.0, 0.0 ));
    joinInGrid(particles[particles.size() - 1]);
    particles.push_back(new Particle( 0.4, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0 ));
    joinInGrid(particles[particles.size() - 1]);*/
    
    double length = 0;
    if (scale)
        length = 0.2;
    double vx = 0.05;
    for (double y = 0.2; y < 0.4+ length; y += 0.03) 
    {
        for (double x = 0.4;  x< 0.6+ length; x += 0.03) 
        {
            //cout << "x is: " << x << " y is: " << y << endl;
            Particle* p=new Particle();
            p->x = x;
            p->y = y;
            p->vx = vx;
            p->vy = -0.2;
            p->ax = 0.0;
            p->ay = 0.0;
            p->rho = 0.0;
            p->p = 0.0;
            particles.push_back(p);
            joinInGrid(particles[particles.size()-1]);
        }
    }
    printGrid();
    
}

// Main function
int SPH_2D(int argc, char** argv,bool scale,bool grid) {
    gridOptimize = grid;
    //初始化所有粒子
    initializeParticles(scale);
    //为每一个粒子计算其密度和压力，
    calculateDensityAndPressure();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winSize, winSize);
    glutCreateWindow("SPH Simulation");
    glutDisplayFunc(display);
    initOpenGL();
    glutTimerFunc(1, update, 0);
    glutMainLoop();
    
    return 0;
}