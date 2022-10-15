clear all

% parameters
parms.R1 = 0.015; % inner beam radius (m) 
parms.R2 = 0.02; % outer beam radius (m) 
parms.E = 181*10^9; % young's modulus of CFRP (Pa) 
parms.I = pi/4*(parms.R2^4-parms.R1^4); % annulus 4th inertia moment (m^4) 
parms.rho_infty = 1.225; % air density (kg/m^3) 
parms.u_infty = 15; % air speed (m/s) 
parms.b = 1.15; % tail wing span (m) 
parms.c = 0.35; % tail chord length (m) 
parms.cl_alpha = 7/0.8; % coeff. of lift derivative wrt alpha (no dim.)
% taken from NACA 0012-34 data 
parms.J = pi/2*(parms.R1^4 - parms.R2^4); % torsion stiffness of beam (Nm/rad)
parms.G = 4.12*10^9; % shear modulus of elasticity of CFRP (Pa) 

parms.P = lift_force(deg2rad(2), parms); % point load (N) 
vertical_deflection(3, parms)
tip_slope(3, parms)
lift_force(deg2rad(2), parms)

function [L, R] = tailwing_force(alpha, parms)
% L: Find lift force (N) on whole tail wing, if oriented with an angle of 
% attack of alpha (rad). 
% R: Find rot. force due to tail wing torque (N*m). 
L = 0.5*parms.rho_infty*parms.u_infty^2*parms.cl_alpha*alpha*parms.b*parms.c;
R = T/x; 
end 

function y = vert_defl(x, parms)
% Find vert. defl. y (m) of two cantilever beams of length x (m), 
% each with one end fixed, the other end free and subject to vert. pt.
% load P. Assume beams have annular cross-section.
y = parms.P*x^3/(3*parms.E*parms.I); 
end 

function phi = tip_slope(x, parms)
% Find the angular slope (rad) of the free end of the same cantilever beam 
% of length x (m).
phi = parms.P*x^2/(2*parms.E*parms.I);
end 

function theta = ang_rot(T, x, parms)
% Calculate angular rotation (rad) of tail wing, as a result of the 
% torque T (N*m) exerted on it. 
theta = 0.5*T*x/(parms.G*parms.J);
end 