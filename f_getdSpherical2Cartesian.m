% https://en.wikipedia.org/wiki/N-sphere
function [x, J] = f_getdSpherical2Cartesian(r, phi)

sin_phi = sin(phi);
x = r .* [cos(phi); 1];

% Transformation
for k = 1:length(phi)
    x(k+1) = x(k+1) * prod(sin_phi(1:k));
end

% Jacobian Matrix
p = [r; phi];
J = jacobian(x, p);

end