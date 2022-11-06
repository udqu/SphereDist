function Kh = f_getdSphereSurfaceKernel(dimension)

% smoothing parameter definition
syms h;
assume(h, {'real','positive'});

% cartesian coordinates in spherical variables
syms r;
phi = sym('phi', [dimension-1,1]);
% assumptions on the spherical variables
assume(r, {'real','positive'});
assume(phi, 'real');
% get the coordinates
[x, J] = f_getdSpherical2Cartesian(r, phi); 

% Surface Measures for the d-Spheres
Sd = [                   2;     % d = 1
                  (2*pi*r);     % d = 2
                (4*pi*r^2);     % d = 3
              (2*pi^2*r^3);     % d = 4
            (8/3*pi^2*r^4);     % d = 5
                (pi^3*r^5);     % d = 6
          (16/15*pi^3*r^6);     % d = 7
            (1/3*pi^4*r^7);     % d = 8
         (32/105*pi^4*r^8);     % d = 9
           (1/12*pi^5*r^9);     % d = 10
        (64/945*pi^5*r^10);     % d = 11
          (1/60*pi^6*r^11);     % d = 12
     (128/10395*pi^6*r^12);     % d = 13
         (1/360*pi^7*r^13);     % d = 14
    (256/135135*pi^7*r^14);     % d = 15
];      

% Integration Values
% d = 2                (2*pi*r) * exp(-r^2/h^2)
% d = 3              (4*pi*r^2) * exp(-r^2/h^2)
% d = 4            (2*pi^2*r^3) * exp(-r^2/h^2)
% d = 5          (8/3*pi^2*r^4) * exp(-r^2/h^2)
% d = 6              (pi^3*r^5) * exp(-r^2/h^2)
% d = 7        (16/15*pi^3*r^6) * exp(-r^2/h^2)
% d = 8          (1/3*pi^4*r^7) * exp(-r^2/h^2)
% d = 9       (32/105*pi^4*r^8) * exp(-r^2/h^2)
% d = 10        (1/12*pi^5*r^9) * exp(-r^2/h^2)
% d = 11     (64/945*pi^5*r^10) * exp(-r^2/h^2)
% d = 12       (1/60*pi^6*r^11) * exp(-r^2/h^2)
% d = 13  (128/10395*pi^6*r^12) * exp(-r^2/h^2)
% d = 14      (1/360*pi^7*r^13) * exp(-r^2/h^2)
% d = 15 (256/135135*pi^7*r^14) * exp(-r^2/h^2)
fprintf('Using Gaussian Surface Kernel...\n');
K_h = (1/Sd(dimension)) * exp(-(x.'*x - r^2)/(h^2));


% create a struct for the kernel
Kh.ndim        = dimension;
Kh.varx        = x;
Kh.varh        = h;
Kh.varp        = phi;
Kh.s_Kh        = K_h;
Kh.s_J         = J;
Kh.s_detJ      = simplify( det(J) );
Kh.f_Kh        = matlabFunction(K_h);
Kh.f_grad_Kh_phi = matlabFunction(simplify( gradient(K_h, phi) ));
Kh.f_diff_Kh_phi = matlabFunction(simplify( diff(K_h, h) ));
% K_h.f_Lx_Kh     = matlabFunction(Lx_Kh);
% K_h.f_Lh_Kh     = matlabFunction(Lh_Kh);

end