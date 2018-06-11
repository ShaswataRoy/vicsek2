/*
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
double eta=0.1;
double L=31;
int N=4000;
double lx=0.2,ly=0.1;

//double ordall[MAX_TIME];
FILE *fp;


//Initial State
double x[MAX_SIZE],y[MAX_SIZE],
    theta[MAX_SIZE],theta_new[MAX_SIZE];
double order,order_new;

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
    order_new=1;
}

//Functions

double old_dist(int i,int j)
{
    return(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
}

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

int permissible(int i,int j)
{
    double l = dist(i,j)/2;
    double A = (l/ly)*(l/ly)-1;
    double det = 4*lx*lx-4*l*l*A;
    if(det>=0)
    {
        double a = (2*lx-sqrt(det))/(2*A);
        double alpha = atan(ly/a);
        if(a>0)
        {
            if(fabs(theta[i]-theta[j])>alpha)
                return 0;
        }
    }
    return 1;
}

void update()
{
    int k;
    double theta_near;
    for(int i=0;i<N;i++)
    {
        theta_near = 0;
        k=0;

        for(int j=0;j<N;j++)
        {
            if(dist(i,j)<1.0)
            {
                theta_near+=theta[j];
                k+=1;
            }
        }
        if(k!=0)
        {
            //theta_new[i] = theta_near/k + eta*uniform(-M_PI,M_PI);
            theta_new[i] = theta_near/k + uniform(-eta/2,eta/2);
        }
    }

    for(int i=0;i<N;i++)
    {
        theta[i]=theta_new[i];
        x[i]=fmod((x[i]+v*cos(theta[i])),L);
        y[i]=fmod((y[i]+v*sin(theta[i])),L);
    }
}

double array_mean(double *array,int min,int max)
{
    double sum=0.0;
    for(int i=min;i<max;i++)
    {
        sum+=array[i];
    }
    return(sum/(max-min));
}

double array_std(double *array,int min,int max)
{
    double sum=0.0;
    double mean=array_mean(array,min,max);
    for(int i=min;i<max;i++)
    {
        sum+=(array[i]-mean)*(array[i]-mean);
    }
    return(sqrt(sum/(max-min)));
}


double array_binder(double *array,int min,int max)
{
    double second = 0.0;
    double fourth = 0.0;
    for(int i=min;i<max;i++)
    {
        second += pow(array[i],2);
        fourth += pow(array[i],4);
    }
    fourth = fourth/(max-min);
    second = second/(max-min);
    return(1-fourth/(second*second*3));
}

//Plot

void compute_order(int no)
{
    initialize();
    int a,b;
    double eps=0.;
    int repeat =1;
    const int no_points=20;
    const double eta_max=2*M_PI;

    a=200;b=10;eps=0.01;repeat=20;

    double order_arr[no_points],eta_arr[no_points];
    double theta_cos,theta_sin;
    double ordall[MAX_TIME];
    int t=0;
    for(int i=0;i<no_points;i++)
    {
        eta_arr[i] = eta_max*i/(1.0*no_points);
    }

    fprintf(fp,"#eta\torder\n");

    for(int r=0;r<repeat;r++)
    {
        printProgress(1.0*r/repeat);


        for(int i=0;i<no_points;i++)
        {
            initialize();
            eta=eta_arr[i];
            for(t=0;t<MAX_TIME;t++)
            {
                if(t>a)
                {
                    if(fabs(array_mean(ordall,t-2*b,t-b)-
                            array_mean(ordall,t-b,t))<eps)
                        break;
                }
                theta_cos = 0.;
                theta_sin = 0.;
                update();

                for(int j=0;j<no;j++)
                {
                    theta_cos+=cos(theta[j]);
                    theta_sin+=sin(theta[j]);
                }
                order_new=sqrt(theta_cos*theta_cos+theta_sin*theta_sin)/no;
                ordall[t]=order_new;
            }

            order_arr[i] += array_mean(ordall,t*3/4,t);
        }
    }
    for(int i=0;i<no_points;i++)
    {
        fprintf(fp,"%lf\t%lf\n",eta_arr[i],order_arr[i]/repeat);
    }

}

int main(int argc,char* argv[])
{

    struct timeval start, end;
    gettimeofday(&start, NULL);
    cout<<"\nN = "<<argv[1]<<"\n";
    N=atoi(argv[1]);L=atof(argv[2]);

    // benchmark code
    string str = "vicsek";
    str+=to_string(N)+".txt";

    fp = fopen(str.c_str(),"w");

    compute_order(N);
    fclose(fp);
    gettimeofday(&end, NULL);

    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;
    cout<<"\rElasped time is "<< delta<<" seconds.\n";
    return(0);
}
*/
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
const int MAX_TIME=500;
double v=0.3;
double eta=0.1;
double L;int N;
double lx=0.9,ly=0.3;

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

double dist(int i,int j)
{
    return(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
}

int permissible(int i,int j)
{
    double l = dist(i,j)/2;
    double A = (l/ly)*(l/ly)-1;
    double det = 4*lx*lx-4*l*l*A;
    if(det>=0)
    {
        double a = (2*lx-sqrt(det))/(2*A);
        double alpha = atan(ly/a);
        if(a>0)
        {
            if(fabs(theta[i]-theta[j])>alpha)
                return 0;
        }
    }
    return 1;
}


void update()
{
    int k;
    double theta_near;
    for(int i=0;i<N;i++)
    {
        theta_near = 0;
        k=0;

        for(int j=0;j<N;j++)
        {
            if(dist(i,j)<1.0 && permissible(i,j)==1)
            {
                theta_near+=theta[j];
                k+=1;
            }
        }
        if(k!=0)
        {
            //theta_new[i] = theta_near/k + eta*uniform(-M_PI,M_PI);
            theta_new[i] = theta_near/k + uniform(-eta/2,eta/2);
        }
    }

    for(int i=0;i<N;i++)
    {
        theta[i]=theta_new[i];
        x[i]=fmod((x[i]+v*cos(theta[i])),L);
        y[i]=fmod((y[i]+v*sin(theta[i])),L);
    }
}

double array_mean(double *array,int min,int max)
{
    double sum=0.0;
    for(int i=min;i<max;i++)
    {
        sum+=array[i];
    }
    return(sum/(max-min));
}

double array_std(double *array,int min,int max)
{
    double sum=0.0;
    double mean=array_mean(array,min,max);
    for(int i=min;i<max;i++)
    {
        sum+=(array[i]-mean)*(array[i]-mean);
    }
    return(sqrt(sum/(max-min)));
}


//Plot

void compute_order(int no)
{
    initialize();
    int a,b;
    double eps=0.;
    int repeat =1;
    const int no_points=20;
    const double eta_max=2*M_PI;

    a=200;b=10;eps=0.01;repeat=200;

    double order_arr[no_points],eta_arr[no_points];
    double theta_cos,theta_sin;
    double ordall[MAX_TIME];
    int t=0;
    for(int i=0;i<no_points;i++)
    {
        eta_arr[i] = eta_max*i/(1.0*no_points);
    }

    fprintf(fp,"#eta\torder\n");

    for(int r=0;r<repeat;r++)
    {
        printProgress(1.0*r/repeat);


        for(int i=0;i<no_points;i++)
        {
            initialize();
            eta=eta_arr[i];
            for(t=0;t<MAX_TIME;t++)
            {
                if(t>a)
                {
                    if(fabs(array_mean(ordall,t-2*b,t-b)-
                            array_mean(ordall,t-b,t))<eps)
                        break;
                }
                theta_cos = 0.;
                theta_sin = 0.;
                update();

                for(int j=0;j<no;j++)
                {
                    theta_cos+=cos(theta[j]);
                    theta_sin+=sin(theta[j]);
                }
                order=sqrt(theta_cos*theta_cos+theta_sin*theta_sin)/no;
                ordall[t]=order;
            }

            order_arr[i] += array_mean(ordall,t*3/4,t);
        }
    }
    for(int i=0;i<no_points;i++)
    {
        fprintf(fp,"%lf\t%lf\n",eta_arr[i],order_arr[i]/repeat);
    }

}

int main()
{

    struct timeval start, end;
    gettimeofday(&start, NULL);
    N=40;L=3.14;

    // benchmark code
    string str = "vicsek";
    str+=to_string(N)+"_09_01.txt";

    fp = fopen(str.c_str(),"w");

    compute_order(N);
    fclose(fp);
    gettimeofday(&end, NULL);

    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;
    cout<<"Elasped time is "<< delta<<" seconds.\n";
    return(0);
}
