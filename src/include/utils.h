#pragma once
#ifndef UTIL_H
#define UTIL_H

#include<vector>
#include<iostream>


//winSizeΪ600����ЧΪ6cm
//�⻬�˰뾶��0.0015��
static double PI = 3.141;
static double winSize = 600;
static double H = 0.0015;
static double MASS = 0.0001;
//��׼�ܶ�
static double RHO0 = 1;
static double K = 20;
static double MU = 0.1;
static double G = -9.8;

static double Vel_atten = 0.2;

//ʱ�䲽
static double timeStep = 0.0001;
static double timeElapsed = 0.0;

const double DEG_TO_RAD = 3.141592653589793 / 180.0;

struct Particle {
    //���ӵ�λ�ã��ٶȣ����ٶ�
    double x, y;
    double vx, vy;
    double ax, ay;
    //�ܶȺ�ѹ��
    double rho;
    double p;
    //������������Ż�
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
    //���ӵ�λ�ã��ٶȣ����ٶ�
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    //�ܶȺ�ѹ��
    double rho;
    double p;
};


//�����񻮷��У���Ϊÿ��Ͱ������ͷ��
struct ListHead
{
    //nextָ���Ͱ�еĵ�һ�����ӣ�tailָ���Ͱ�е����һ������
    Particle* next = nullptr;
    Particle* tail = nullptr;
    int num=0;

    ListHead();
    void addNew(Particle* par);
};

//ȫ������
extern std::vector<std::vector<ListHead*>> grid;

void initGrid(int num);
void joinInGrid(Particle* par);
void printGrid();
void removeFromOld(Particle* par);
void checkFunc();
#endif // UTIL_H

