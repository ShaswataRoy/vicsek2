//To complile this code, type following command.
// gcc vicsek.c -o vicsek -lm && ./vicsek
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "random.h"

using namespace std;

const long long int N=2048;
int steps=300000, near_count;
double x[N], y[N], vx[N], vy[N], theta[N], mean_ang[N], list_ang[N];
double L=64, vc=0.5, dt=1, eta=0.15, r=1;

//double unirand(double a, double b); // fn to return uniform random number in (a,b)
//double dist(double x1, double x2, double y1, double y2, double L); // fn to compute the distance between i th particle and j th particle at periodic boundary condition
//void pdist(double D[], double x[], double y[], int n, double L); // fn to compute the pairwise distance a periodic boundary condition
//void cal_mean_angle(double D[], double mean_ang[], double theta[], int N, double r); // fn to compute mea angle
//void filewrite(double x[], double y[], double vx[], double vy[], int iter, int N); // fn to write down the positions and the velocities

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

double bdist(int i,int j){
    return(dist(x[i], x[j], y[i], y[j]));
}

void cal_mean_angle(){

    int count;

    for (int i=0; i<N; i++){
        int near_count = 0;
        for (int j=0; j<N; j++){
            if (bdist(i,j)<= r){
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

void filewrite(int iter){
        FILE *fp;
        char filename[20];
        sprintf(filename, "./data3/%d.txt", iter);

        fp = fopen(filename, "w");
        if (fp == NULL) {
            printf("cannot open\n");
            exit(1);
        }

        for (int i=0; i<N; i++){
            fprintf(fp, "%lf\n", theta[i]);
        }
        fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int main(void){

    clock_t begin = clock();
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
     double mean_sin, mean_cos;


	// make "data" directory to store positions and velocities data.
	mkdir("data", 0755);

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
    for (int iter=0; iter<steps; iter++){
        if(iter>steps/2)
        {
            filewrite(iter);
        }

        // calculate mean angle
       cal_mean_angle();


        // update theta
        for (int i=0; i<N; i++){
            theta[i] = mean_ang[i] + M_PI * eta * unirand(-0.5, 0.5);
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
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout<<"Time: "<<elapsed_secs<<"\n";
    return 0;
}
