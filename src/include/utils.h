#pragma once

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

struct Particle {
    //���ӵ�λ�ã��ٶȣ����ٶ�
    double x, y;
    double vx, vy;
    double ax, ay;
    //�ܶȺ�ѹ��
    double rho;
    double p;
};
