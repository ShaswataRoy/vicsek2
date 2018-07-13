// To complile this code, type following command.
// gcc vicsek.c -o vicsek -lm && ./vicsek
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "progress.h"
#include "random.h"

using namespace std;

FILE *fp;

double unirand(double a, double b); // fn to return uniform random number in (a,b)
const int N=100;
double L=5;
double r =1.0;
double order = 0.;
double eta = 3.;
double vc=0.1, dt=1;
double x[N],y[N],vx[N],vy[N],D[N*N],theta[N], mean_ang[N], list_ang[N];
//double dist(); // fn to compute the distance between i th particle and j th particle at periodic boundary condition
void pdist(); // fn to compute the pairwise distance a periodic boundary condition
void cal_mean_angle(); // fn to compute mea angle
void filewrite(int iter); // fn to write down the positions and the velocities

double unirand(double a, double b){
    return a + (b-a) * Uniform();
}

double dist(double x1, double x2, double y1, double y2){
    double dx, dy;
    dx = fabs(x1 - x2);
    if (dx > L-dx)
        dx = L-dx;
    dy = fabs(y1 - y2);
    if (dy > L-dy)
        dy = L-dy;
    return sqrt(dx*dx + dy*dy);
}

void pdist(){
   int i, j, count=0;
   for (i = 0; i < N; i++){
      for (j = 0; j < N; j++){
         D[count] = dist(x[i], x[j], y[i], y[j]);
         count += 1;
      }
   }
}

void cal_mean_angle(){
    double A[N][N], list_ang[N];
    int count;

    // D[N*N] --> A[N][N]
    count = 0;
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            A[i][j] = D[count];
            count += 1;
        }
    }

    for (int i=0; i<N; i++){
        int near_count = 0;
        for (int j=0; j<N; j++){
            if (A[i][j] <= r){
                list_ang[near_count] = theta[j];
                near_count += 1;
            }
        }
        double mean_sin=0.0, mean_cos=0.0;
        for (int j=0; j<near_count; j++){
            if(cos(theta[i]-list_ang[j])>0)
            {
                mean_sin += sin(list_ang[j]);
                mean_cos += cos(list_ang[j]);
            }
            else
            {
                mean_sin -= sin(list_ang[j]);
                mean_cos -= cos(list_ang[j]);
            }
        }

        mean_sin = mean_sin / near_count;
        mean_cos = mean_cos / near_count;

        mean_ang[i] = atan2(mean_sin, mean_cos);
    }
}

double compute_order()
{
    // pairwise distance
    pdist();
     // calculate mean angle
    cal_mean_angle();


     // update theta
     for (int i=0; i<N; i++){
         theta[i] = mean_ang[i] + eta * unirand(-0.5, 0.5);
     }

     // update positions and velocities
     for (int i=0; i<N; i++){
         vx[i] = vc * cos(theta[i]); vy[i] = vc * sin(theta[i]);
         x[i] += vx[i] * dt; y[i] += vy[i] * dt;
     }

     // periodic boundary condition
     for (int i=0; i<N; i++){
         x[i] = fmodl(x[i] + L, L);
         y[i] = fmodl(y[i] + L, L);
     }

    double sum_cos = 0.0;
    double sum_sin = 0.0;

    for (int i=0; i<N; i++){
        sum_cos += cos(2*theta[i]);
        sum_sin += sin(2*theta[i]);
    }

    order = sqrt(sum_cos*sum_cos+sum_sin*sum_sin)/N;

    return(order);
}

void filewrite(){
        fprintf(fp, "%.4lf\t%.4lf\n",eta, order);
}

void generate()
{
    int steps=500, near_count;
    double mean_sin, mean_cos;
    double eps = 0.001;
    // seed
    init_genrand((unsigned)time(NULL));

    for(int i=0; i<N; i++){
        theta[i] = unirand(-M_PI, M_PI);
        x[i] = unirand(0.0, L);
        y[i] = unirand(0.0, L);
        vx[i] = vc * cos(theta[i]);
        vy[i] = vc * sin(theta[i]);
    }

    //
    // start vicsek model
    //
    int iter=0;

    while(iter<steps){
        iter++;
        compute_order();
    }

    double mean_order = 0.;

    for(int i=0;i<100;i++)
    {
        mean_order += compute_order();
    }

    order = mean_order/100;

    filewrite();
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int main(void){
	/*
	 * N: the number of particle
	 * steps: the number of run steps
	 * x[i], y[i]: the x,y-position of i th particle
	 * vx[i], vy[i]: the vx,vy-velosity of i th particle
	 * theta[i]: the angle of i th particle
	 * D[i][j]: pairwise distance between i th particle and j th particle
	 *
	 * L: system size
	 * vc: speeds of the particles
	 * dt: time interval
	 * eta: order
	 * r: interaction radius
	 *
	 */

     string str = "vicsek";
     str+=to_string(N)+"n.txt";

     fp = fopen(str.c_str(),"w");

     for(eta=0.2;eta<=6.0;eta+=0.2)
     {
         printProgress(eta/6.0);
         generate();
     }

     fclose(fp);

    return 0;
}
