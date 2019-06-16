
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <random>
#include <cstdlib>
#include <iostream>
#include "gnuplot.h"

using namespace std;
using namespace MPI;

double func(double x);
void NewtonDiffTable(int npts, double *xpts, double *funcvals,double * newton_coeffs);
void   CreateGrid_EvenlySpaced(int npts, double *x,double a, double b);
double NewtonInterpolant(double x, int npts, double * xpts,double * newton_coeffs);
double NewtonDiffFunction(int start_index, int ending_index,double * xpts, double * funcvals);



int main(int argc, char * argv[]){
    int i;
    int degree, polypnts;
    int npts = 40000;
    int my_rank, comm_size;
    
     ofstream myfile;
     myfile.open("out.txt");
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,& my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    
    //number of points used for plotting
    double xpt, soln, approx;
    cout << "Enter the degree of the interpolating polynomial: ";
    cin >> degree;
    polypnts = degree+1; //number of points is equal to 1 + degree
    
    double * poly_xpts = new double[polypnts];
    double * func_vals = new double[polypnts];
    double * newton_coeffs = new double[polypnts];
     double* x =new double[npts]();
    CreateGrid_EvenlySpaced(polypnts, poly_xpts, -1.0, 1.0);
    
    int ppn= npts/comm_size;
    
    if(my_rank !=0){
         
        for(i=0;i<polypnts;i++){        
            func_vals[i] = func(poly_xpts[i]);            
        }
        MPI_Send(func_vals, ppn,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        NewtonDiffTable(polypnts, poly_xpts, func_vals,newton_coeffs);  
    }
    else{
        for(i=0;i<ppn;i++){ 
            xpt = -1.0 + 2.0*i/(npts-1);       
            x[i]=func(xpt);
            for(int j=1;j<comm_size;j=j+1){
                MPI_Recv(x,ppn,MPI_DOUBLE,MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            approx = NewtonInterpolant(xpt, polypnts,poly_xpts, newton_coeffs);
            cout << xpt << " " << soln << " " << approx << endl;
            myfile<<xpt << " " << soln << " " << approx <<"\n";   
        }
    }

     //plots
   // GnuplotPipe gp;
    //gp.sendLine("plot sin(x)  lw 4");
    
    
    delete[] x;
    delete[] poly_xpts;
    delete[] func_vals;
    delete[] newton_coeffs;
    
    myfile.close();
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}





double func(double x){
    double y;
    y = 1.0 + 25.0*x*x;
    y = 1.0/y;
    return y;
}


double NewtonInterpolant(double x, int npts, double * xpts,double * newton_coeffs){
    int i,j;
    double sum = 0.0, xval;
    for(i=0;i<npts;i++){
        xval = 1.0;
        for(j=0;j<i;j++)
            xval = xval*(x-xpts[j]);
        sum = sum + newton_coeffs[i]*xval;
    }
    return sum;
}


void NewtonDiffTable(int npts, double *xpts, double *funcvals,double * newton_coeffs){
    int i,j;
    for(i=0;i<npts;i++)
        newton_coeffs[i] = NewtonDiffFunction(0,i, xpts, funcvals);
}


double NewtonDiffFunction(int start_index, int ending_index,double * xpts, double * funcvals){
    double val;
    int diff = ending_index-start_index;
    if(diff == 0){
        val = funcvals[start_index];
    }
    else{
    val = (NewtonDiffFunction(start_index,ending_index-1,xpts,funcvals) -
            NewtonDiffFunction(start_index+1,ending_index,xpts,funcvals))/(xpts[start_index]-xpts[ending_index]);
}
return val;
}

void CreateGrid_EvenlySpaced(int npts, double *x, double a, double b){
    double dx = (b-a)/(npts-1.0);
    for(int i=0;i<npts;i++)
        x[i] = a + i*dx;
    return;
}







