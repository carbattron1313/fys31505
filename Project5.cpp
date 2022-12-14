//
// Project5.hpp
//
//
// Created by Alejandro Carballido Mantecón, Elena Muñoz Rivas, David Martínez Hernández and Antonio Gómez Garrido.
//
// Program that creates a simulate the two-dimensional time-dependent Schrödinger equation, and use it to study a double-slit-in-a-box setup.



#include "schro.hpp"


    
int main(){
    
    //Define all the parameters for the simulation
    double h = 0.005;
    double dt = 2.5e-5;
    double T = 0.002;
    double xc = 0.25;
    double yc = 0.5;
    double sigma_x = 0.05;
    double sigma_y = 0.2;
    double px = 200;
    double py = 0;
    double v0 = 1e+10;
    int M = 1./h;
    int n_slits = 3;


    //Set up the potential
    arma::sp_mat potential = V(v0, h, M, n_slits);
    
    //Define the matrices A and B
    arma::sp_cx_mat A = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));
    arma::sp_cx_mat B = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));


    //Fill matrices A and B
    mat_fill(M,h,dt,potential,A,B);

    //Compute the initial state of u_0
    arma::cx_mat u_0 = set_up_u_0(xc, yc, sigma_x, sigma_y, px, py, M, h);

    //Define a cube for storing all data for different times
    arma::cx_cube time_cube = arma::cx_cube(M-2,M-2,T/dt);
    time_cube.slice(0) = u_0;

    //Compute the wave function for each time and add it to the data cube
    for (int i=1; i<T/dt; i++){
        arma::cx_mat u_1 = time_cube.slice(i-1);
        arma::cx_vec u_2 =  u_2(A,B,u_1.as_col());
        arma::cx_mat u_2_mat = reshape(u_2,M-2,M-2);
        time_cube.slice(i) = u_2_mat;
    }
    
    //save data cube
    time_cube.save("cube2_3s.bin");
    return 0;
    
}
