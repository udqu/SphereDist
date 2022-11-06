function integral_f = f_getMultiDimensionalSphereSurfaceIntegral(Kh)

% initialize the integral with Jacobian times the function
integral_f = abs(Kh.s_detJ) * Kh.s_Kh;

% perform the integration
for k = 1:length(Kh.varp)-1
    integral_f = int(integral_f, Kh.varp(k), -pi/2, pi/2);
end
integral_f = simplify(expand( int(integral_f, Kh.varp(end), -pi, pi) ));

end