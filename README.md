Time dependent Schrödinger equation for a double slit in a box
------------------------
There are one .cpp file, Project5.cpp, and one .hpp file, Project5.hpp.

Project5.hpp
--------------
Includes the functions:
- k: converts the indices of the matrix into a vector index.
- mat_construct: constructs the matrices A and B following the Crank-Nicolson method
- mat_fill: fills the matrices diagonals and subdiagonals
- u_next: obtains the wave function in the next time step
- set_up_u_0: sets up the initial wave function
- V: sets up the slit barrier 



Project5.cpp
--------------
Main program to to simulate the two-dimensional time-dependent Schrödinger equation, and use it to study a double-slit-in-a-box setup.

Build: g++ -O3 Project5.cpp -std=c++14 -larmadillo -o Project5.exe
Run: ./Project5.exe

Schrodinger_5.py
--------------
Python script that reads the datas obtained in the C++ and generates plots of the the analytical results and the values obtained. It contains the code of the different plots made. The desviation of the total probability from 1.0 as a function of time with a log plot in the y axis. The three colourmap plots that illustrate the time evolution of the 2D probability function for different times with the real and imaginary parts. The detection probability along the particle with a detector screen x=0.8 and time t=0.002. The one-dimensional probability functionis normalise. 

To run the program is necessary to install pyarma.

Run command: python3 Isingmodel.py
