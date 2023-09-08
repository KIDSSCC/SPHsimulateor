#pragma once
#ifndef UTIL_H
#define UTIL_H

#include<vector>
#include<iostream>


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

const double DEG_TO_RAD = 3.141592653589793 / 180.0;

struct Particle {
    //粒子的位置，速度，加速度
    double x, y;
    double vx, vy;
    double ax, ay;
    //密度和压力
    double rho;
    double p;
    //用于网格采样优化
    int i = -1, j = -1;
    Particle* next = nullptr;
    Particle* last = nullptr;
    Particle(double x, double y, double vx, double vy, double ax, double ay, double rho, double p)
        : x(x), y(y), vx(vx), vy(vy), ax(ax), ay(ay), rho(rho), p(p)
    {
        
    }
    Particle() {};

};



struct Particle3D {
    //粒子的位置，速度，加速度
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    //密度和压力
    double rho;
    double p;
};


//在网格划分中，作为每个桶的链表头；
struct ListHead
{
    //next指向该桶中的第一个粒子，tail指向该桶中的最后一个粒子
    Particle* next = nullptr;
    Particle* tail = nullptr;
    int num=0;

    ListHead();
    void addNew(Particle* par);
};

//全局网格
extern std::vector<std::vector<ListHead*>> grid;

void initGrid(int num);
void joinInGrid(Particle* par);
void printGrid();
void removeFromOld(Particle* par);
void checkFunc();
#endif // UTIL_H

