clc; clear; close all;

% Set the parameters for the problem
d = 3;     % dimension of the problem
N = 9;     % number of points on d-dimensional sphere
R = 100;   % radius of the d-dimensional sphere

% Get the d-Sphere Surface Kernel (default Gaussian)
Kh = f_getdSphereSurfaceKernel(d);

