//
// Project5.hpp
//  
//
// Created by Alejandro Carballido Mantecón, Elena Muñoz Rivas, David Martínez Hernández and Antonio Gómez Garrido.
//

#ifndef Schrodinger_Project5_hpp
#define Schrodinger_Project5_hpp

#include <stdio.h>
#include <iostream>
#include <complex>
#include <armadillo>

//Needed to use for the complex numbers
using namespace std::complex_literals;


//Function that calculates the index k with the indices i and j
inline int k(int i, int j, int M){
    return i + (j-1)*(M-2) -1;
}


//Function that constructs a matrix with values according to Crank-Nicolson:
void mat_construct(int M, arma::cx_double r, arma::cx_vec main_diag, arma::sp_cx_mat& matrix){
    
    matrix.diag(0) = main_diag;
    matrix.diag(1).fill(r);
    matrix.diag(-1).fill(r);
    matrix.diag(M-2).fill(r);
    matrix.diag(-(M-2)).fill(r);

    for (int i=1; i<= matrix.diag(1).n_rows; i++){
        if (i%(M-2) == 0){
            matrix.diag(1)(i-1) = 0.;
            matrix.diag(-1)(i-1) = 0.;
        }
    }
}


//Function that fill both A and B according to Crank-Nicolson:
void mat_fill(int M, double h, double dt, arma::sp_mat V, arma::sp_cx_mat& mat_A, arma::sp_cx_mat& mat_B){
    arma::cx_vec a = arma::cx_vec((M-2)*(M-2));
    arma::cx_vec b = arma::cx_vec((M-2)*(M-2));
    arma::cx_double r = 1i*dt/(2.*h*h);

    for (int i=1; i < (M-2); i++){
        for (int j=1; j < (M-2); j++){
            a(k(i,j,M)) = 1. + 4.*r + 1i*dt*((double) V(i,j))/2.;
            b(k(i,j,M)) = 1. - 4.*r - 1i*dt*((double) V(i,j))/2.;
        }
    }
    mat_construct(M,-r,a, mat_A);
    mat_construct(M,r,b, mat_B);
}


//Function that obtains the wave function at n+1 using the armadillo solver
arma::cx_vec u_next(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec u){
    arma::cx_vec b = B*u;
    arma::cx_vec u_next = arma::spsolve(A,b);
    return u_next;
}

//Function that sets up the initial state u0
arma::cx_mat set_up_u_0(double xc, double yc, double sigma_x, double sigma_y, double px, double py, double M, double h){
    arma::cx_mat u_0 = arma::cx_mat(M,M);
    //compute the denominator here to avoid doing it inside the loop
    double subx=2*sigma_x*sigma_x;
    double suby=2*sigma_y*sigma_y;
    for (int i=1; i<M-2; i++){
        for (int j=1; j<M-2; j++){
            u_0(j,i) = exp(-(i*h-xc)*(i*h-xc)/(subx) - (j*h-yc)*(j*h-yc)/(suby) + 1i*px*(i*h-xc) + 1i*py*(j*h-yc));
        }
    }
    
    //Remove the boundaries of u_0
    arma::uvec i_row = {0,199};
    arma::uvec i_col = {0,199};
    u_0.shed_cols(i_row);
    u_0.shed_rows(i_col);
    
    return u_0/std::sqrt(arma::sum(u_0.as_col()%arma::conj(u_0.as_col())).real());
}


//Function that initialized the potential V and constructs the barrier with all the slits you want.
arma::sp_mat V(double v0, double h, double M, double n_slits){
    int wall_thickness = 0.02/h;
    int wall_position = 0.5/h-0.5*wall_thickness; //set the start half barrier to the left
    int slit_aperture = 0.05/h;
    int slit_separation = 0.05/h;
    //check that you didn´t ask for too many slits or for a bigger separation than the size of the box
        
    int slits_total_width = slit_aperture*n_slits+slit_separation*(n_slits-1);
        
    if (slits_total_width*h>1){
        std::cout<<"Try to reduce the slit separation or the number of slits"<<std::endl;
        std::exit(1);
    }
        
    int slits_begin = 0.5/h-0.5*slits_total_width;
        
    arma::sp_mat V = arma::sp_mat((M-2),(M-2)); //we don't consider the boundaries

    for (int i=0; i<M-2; i++){
        for (int j=0; j<M-2; j++){
            //set the wall in the x direction
            if (wall_position <= i && i <= wall_position+wall_thickness){
                V(j,i) = v0;
                if (slits_begin<=j && j<=(slits_begin+slits_total_width)){ //This is just for avoiding looping when we are not at the slits
                    for (int n=0;n<n_slits;n++){
                        if (slits_begin+n*(slit_aperture+slit_separation) <= j && j <= slits_begin+(n+1)*slit_aperture+n*slit_separation){
                                V(j,i) = 0;
                        }
                    }
                }
            }
        }
    }
    return V;
}
#endif /* schro_hpp */
