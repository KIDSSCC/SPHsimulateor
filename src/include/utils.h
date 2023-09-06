#pragma once

//winSize为600，等效为6cm
//光滑核半径，0.0015米
static double PI = 3.141;
static double winSize = 600;
static double H = 0.0015;
static double MASS = 0.0001;
//标准密度
static double RHO0 = 1;
static double K = 20;
static double MU = 0.1;
static double G = -9.8;

static double Vel_atten = 0.2;

//时间步
static double timeStep = 0.0001;
static double timeElapsed = 0.0;

struct Particle {
    //粒子的位置，速度，加速度
    double x, y;
    double vx, vy;
    double ax, ay;
    //密度和压力
    double rho;
    double p;
};
