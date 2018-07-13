#include<stdio.h>
#include <cstring>
#include<math.h>
#include <sys/time.h>
#include <random>
#include <iostream>
#include <omp.h>
#include "progress.h"


using namespace std;

//Define Variables
const int MAX_SIZE=10000;
const int MAX_TIME=10000;
double v=0.03;
double eta=2.0;
double L=10;
int N=40;
double r = 1.0;
int num_steps;

//double ordall[MAX_TIME];
FILE *fp;


//Initial State
double x[MAX_SIZE],y[MAX_SIZE],
    theta[MAX_SIZE],theta_new[MAX_SIZE];
double order;

//Initialize
double uniform(double a,double b)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(a, b);
    return(dis(gen));
}

void initialize_array(double *array,double a,double b)
{
    for(int i=0;i<N;i++)
    {
        array[i] = uniform(a,b);
    }
}

void initialize()
{
    initialize_array(x,0,L);
    initialize_array(y,0,L);
    initialize_array(theta,-M_PI,M_PI);
    initialize_array(theta_new,-M_PI,M_PI);

    order=0;
}

//Functions
double dist(int i,int j){
    double dx, dy;
    dx = fabs(x[i] - x[j]);
    if (dx > L-dx)
        dx = L-dx;
    dy = fabs(y[i] - y[j]);
    if (dy > L-dy)
        dy = L-dy;
    return sqrt(dx*dx + dy*dy);
}

void update()
{
    double theta_near_cos,theta_near_sin;
    for(int i=0;i<N;i++)
    {
        theta_near_cos = 0;
        theta_near_sin = 0;

        for(int j=0;j<N;j++)
        {
            if(dist(i,j)<r)
            {
                if(cos(theta[i]-theta[j])>0)
                {
                    theta_near_cos += cos(theta[j]);
                    theta_near_sin += sin(theta[j]);
                }
                else
                {
                    theta_near_cos -= cos(theta[j]);
                    theta_near_sin -= sin(theta[j]);
                }
            }
        }
            //theta_new[i] = theta_near/k + eta*uniform(-M_PI,M_PI);
        theta_new[i] = atan2(theta_near_sin,theta_near_cos) +
                            uniform(-eta/2,eta/2);
    }

    for(int i=0;i<N;i++)
    {
        theta[i]=theta_new[i];
        x[i]=fmod((x[i]+v*cos(theta[i])),L);
        y[i]=fmod((y[i]+v*sin(theta[i])),L);
    }
}

void compute_order()
{
    double sum_cos=0.0;
    double sum_sin=0.0;
    num_steps = 100;
    fprintf(fp,"#eta\torder\n");

    for(eta=0;eta<=2.0;eta+=0.1)
    {
        initialize();
        for(int i=1;i<=num_steps;i++)
        {
            update();
        }

        for(int i=0;i<N;i++)
        {
            sum_cos += cos(2*theta[i]);
            sum_sin += sin(2*theta[i]);
        }
        order = sqrt(sum_cos*sum_cos+sum_sin*sum_sin)/N;
        fprintf(fp,"%lf\t%lf\n",eta,order);
    }

}

int main()
{
    string str = "vicsek";
    str+=to_string(N)+".txt";

    fp = fopen(str.c_str(),"w");
    compute_order();

    return(0);
}
