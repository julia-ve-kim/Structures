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

parms.L = lift_force(deg2rad(2), parms); % point load (N) 
parms.R = rotational_force(10, 3); 
parms.FT1 = total_force_end_1(parms);
parms.FT2 = total_force_end_2(parms); 
parms.y1 = vertical_deflection_end_1(3, parms);
parms.y2 = vertical_deflection_end_2(3, parms);
parms.theta = ang_rotation(10, 3, parms);

total_ang_rotation(parms) 

function [L, R] = tailwing_force(alpha, parms)
% L: Find lift force (N) on whole tail wing, if oriented with an angle of 
% attack of alpha (rad). 
% R: Find rot. force due to tail wing torque (N*m) on either end of beam. 
L = 0.5*parms.rho_infty*parms.u_infty^2*parms.cl_alpha*alpha*parms.b*parms.c;
R = (0.5*T)/(0.5*x);
end 

function [FT1,FT2] = net_force_end(parms) 
% FT1: Find total upwards force (N) exerted on end of beam, where
% rot. and lift forces add constructively (i.e., in same direction). 
% FT2: ... add destructively (i.e., in opposite directions) 
% Suppose lift force exerted on whole wing is experienced equally by 
% both ends.
FT1 = parms.L/2 + parms.R; 
FT2 = parms.L/2 - parms.R; 
end 

function [y1, y2, y_tot] = vert_defl(x, parms)
% y1, y2: Find vert. deflection y (m) of ends 1 and 2, given each is 
% subject to a tot. vert. pt. load FT1 (N). Model is a cantilever beam of 
% length x, with one end fixed, the other end free and subject to a vert. 
% pt. load P. Assume the beam has an annular cross-section.
% y_tot: Find the total vertical deflection of the whole tail wing. 
y1 = parms.FT1*x^3/(3*parms.E*parms.I); 
y2 = parms.FT2*x^3/(3*parms.E*parms.I); 
y_tot = parms.y1 + parms.y2;
end 

function [theta, theta_tot] = ang_rot(T, x, parms)
% theta: Calculate angular rotation (rad) of either end, as a result of the 
% torque T (N*m) exerted on the whole tail wing. 
% theta_tot: Find tot. ang. rot. (rad), experienced by whole tail wing. 
theta = 0.5*T*x/(parms.G*parms.J);
theta_tot = parms.theta*2;
end 