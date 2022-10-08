clear all

% parameters  
parms.R1 = 0.015; % inner beam radius (m) 
parms.R2 = 0.02; % outer beam radius (m) 
parms.E = 14.5 * 10^9; % young's modulus of CFRP (Pa) 
parms.I = pi/4*(parms.R2^4-parms.R1^4); % annulus 4th inertia moment (m^4) 
parms.rho_infty = 1.225; % air density (kg/m^3) 
parms.u_infty = 15; % air speed (m/s) 
parms.b = 1.15; % tail wing span (m) 
parms.c = 0.35; % tail chord length (m) 
parms.cl_alpha = 7/0.8; % coeff. of lift derivative wrt alpha (no dim.)
% taken from NACA 0012-34 data 

parms.P = lift_force(deg2rad(2), parms); % point load (N) 
vertical_deflection(3, parms)
tip_slope(3, parms)
lift_force(deg2rad(2), parms)

function y = vertical_deflection(x, parms)
% Find the vertical deflection y (m) of two cantilever beams of length x (m), 
% each with one end fixed, the other end free and subject to a vertical point
% load P. Assume beams have annular cross-section.

y = 2*parms.P*x^3/(3*parms.E*parms.I); 
end 

function theta = tip_slope(x, parms)
% Find the angular slope (rad) of the free end of the same cantilever beam 
% of length x (m).

theta = parms.P*x^2/(2*parms.E*parms.I);
end 

function L = lift_force(alpha, parms)
% Calculate the lift force (N) on the tail wing, given that it is
% oriented with an angle of attack of alpha (rad). 

L = 0.5*parms.rho_infty*parms.u_infty^2*parms.cl_alpha*alpha*parms.b*parms.c;
end 

