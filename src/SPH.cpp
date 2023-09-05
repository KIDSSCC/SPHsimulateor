#include <GL/glut.h>
#include <cmath>
#include <vector>
#include<iostream>
using namespace std;
// 常量的声明
double PI = 3.141;

//winSize为700，等效为7cm
double winSize = 700;
//光滑核半径，0.02米
double H = 0.02;
double MASS = 1.0;
double RHO0 = 1000.0;
double K = 20;
double MU = 0.1;
double G = -98;

double Vel_atten = 0.7;
double NORMAL_DENSITY = 315 * MASS / (64 * PI * pow(H, 9));
double NORMAL_PRESSURE = -45 * MASS / (PI * pow(H, 6));
double NORMAL_VISCOUS = 45 * MU * MASS / (PI * pow(H, 6));


double kernel_weight(float r, float h) {
    if (r > h) {
        return 0.0f;
    }
    float q = r / h;
    float alpha = 15.0f / (2.0f * PI * h * h);
    float beta = (1.0f - q) * (1.0f - q);
    return alpha * beta * beta * beta;
}
struct Particle {
    //粒子的位置，速度，加速度
    double x, y;
    double vx, vy;
    double ax, ay;
    //密度和压力
    double rho;
    double p;
};

void kernel_gradient(double r, double h, double& gx, double& gy, Particle& p1, Particle& p2) {
    if (r > h) {
        gx = 0.0f;
        gy = 0.0f;
        return;
    }
    float q = r / h;
    float alpha = -45.0f / (PI * h * h * h);
    float beta = (1.0f - q) * (1.0f - q);
    gx = alpha * beta * beta * (p1.x - p2.x) / r;
    gy = alpha * beta * beta * (p1.y - p2.y) / r;
}

std::vector<Particle> particles;
double timeStep = 0.01;
double timeElapsed = 0.0;
double maxX = 1.0, maxY = 1.0;
double minX = 0.0, minY = 0.0;

//以米为单位，计算两个粒子之间的距离
double distance(Particle p1, Particle p2) {
    auto res = sqrt(pow((p1.x - p2.x) * winSize / 100, 2) + pow((p1.y - p2.y) * winSize / 100, 2));
    return res/100;
}

double K_Poly6 = 4 / (PI * pow(H, 8));
double W_Poly6(double r)
{
    double res = K_Poly6 * pow((pow(H, 2) - pow(r, 2)), 3);
    return res;
}

double pressure(Particle p) {
    return 0.4 * (p.rho/RHO0 - 1);
}

// Function to calculate viscosity
double viscosityX(Particle p1, Particle p2) {
    double r = distance(p1, p2);
    return NORMAL_VISCOUS * (p2.vx - p1.vx) / p1.rho * (H - r);
}
double viscosityY(Particle p1, Particle p2) {
    double r = distance(p1, p2);
    return NORMAL_VISCOUS * (p2.vy - p1.vy) / p1.rho * (H - r);
}

// Function to calculate acceleration
void calculateAcceleration(Particle& p) {
    p.ax = 0.0;
    p.ay = G;

    for (Particle &other : particles) {
        if (&p == &other) {
            continue;
        }
        double r = distance(p, other);
        if (r <2*H&&r!=0) {
            
            /*double pressure = p.p;
            double gx = 0.0f;
            double gy = 0.0f;
            kernel_gradient(r, H, gx, gy, p, other);
            p.ax -= MASS * (pressure / pow(p.rho, 2) + other.p / pow(other.rho, 2) + MU * (other.vx - p.vx) / other.rho) * gx;
            p.ay -= MASS * (pressure / pow(p.rho, 2) + other.p / pow(other.rho, 2) + MU * (other.vy - p.vy) / other.rho) * gx;
            */
            
            double pXTerm = NORMAL_PRESSURE * (-(other.x - p.x)*winSize / r * (other.p + p.p) / (2 * other.rho) * pow(H - r, 2));
            double pYTerm= NORMAL_PRESSURE * (-(other.y - p.y) * winSize / r * (other.p + p.p) / (2 * other.rho) * pow(H - r, 2));
            double vXTerm = viscosityX(p, other);
            double vYTerm = viscosityY(p, other);
            p.ax -= (pXTerm + vXTerm)/p.rho;
            p.ay -= (pYTerm + vYTerm)/p.rho;
            
        }
    }
}

// Function to update particle positions and velocities
void updateParticles() {
    int debug = 0;
    for (Particle& p : particles) {
        p.vx += p.ax * timeStep;
        p.vy += p.ay * timeStep;
        p.x += p.vx * timeStep;
        p.y += p.vy * timeStep;
        if (p.x > maxX) {
            p.vx = -Vel_atten *p.vx;
            p.x = maxX;
        }
        if (p.x < minX) {
            p.vx = -Vel_atten *p.vx;
            p.x = minX;
        }
        if (p.y > maxY) {
            p.vy = -Vel_atten *p.vy;
            p.y = maxY;
        }
        if (p.y < minY) {
            p.vy = -Vel_atten *p.vy;
            p.y = minY;
        }
        calculateAcceleration(p);
    }
}

// Function to initialize particles
void initializeParticles() {

    //particles.push_back({ 0.4, 0.6, 6, 0.0, 0.0, 0.0, 0.0, 0.0 });
    //particles.push_back({ 0.4, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    for (double x = 0.1; x < 0.9; x += 0.1) {
        for (double y = 0.1; y < 0.9; y += 0.1) {
            Particle p;
            p.x = x;
            p.y = y;
            p.vx = 6;
            p.vy = -4;
            p.ax = 0.0;
            p.ay = 0.0;
            p.rho = 0.0;
            p.p = 0.0;
            particles.push_back(p);
        }
    }
    
}

// Function to calculate density and pressure for each particle
void calculateDensityAndPressure() {
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i].rho = 0;
        for (int j = 0; j < particles.size(); j++)
        {
            double r = distance(particles[i], particles[j]);
            //考虑与当前粒子两厘米范围内的所有粒子
            if (r < H)
            {
                particles[i].rho += NORMAL_DENSITY * pow((pow(H, 2) - pow(r, 2)), 3);
                //particles[i].rho += MASS * kernel_weight(r, H);
            }
        }
        particles[i].p = pressure(particles[i]);
    }
}

// Function to draw particles
void drawParticles() {
    int numSegment = 200;
    float radius = 0.01;

    glBegin(GL_POINTS);
    for (Particle p : particles) {
        glColor3f(0.0, 0.0, 1.0);
        for (int i = 0; i < numSegment; i++)
        {
            float angle = i * 2.0 * PI / numSegment;
            float x = p.x + radius * cos(angle);
            float y = p.y + radius * sin(angle);
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
    glutTimerFunc(100, update, 0);
}

// Function to initialize OpenGL
void initOpenGL() {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(minX, maxX, minY, maxY);
}

// Main function
int main(int argc, char** argv) {
    //初始化所有粒子
    initializeParticles();
    //为每一个粒子计算其密度和压力，
    calculateDensityAndPressure();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winSize, winSize);
    glutCreateWindow("SPH Simulation");
    glutDisplayFunc(display);
    initOpenGL();
    glutTimerFunc(10, update, 0);
    glutMainLoop();
    
    return 0;
}