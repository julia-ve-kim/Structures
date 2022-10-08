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
parms.J = 12000*pi/180; % torsion stiffness of beam (Nm/rad)
parms.G = 228*10^3; % modulus of elasticity of beam (Pa) 

parms.L = lift_force(deg2rad(2), parms); % point load (N) 
parms.R = rotational_force(10, 3); 
parms.FT1 = total_force_end_1(parms);
parms.FT2 = total_force_end_2(parms); 
parms.y1 = vertical_deflection_end_1(3, parms);
parms.y2 = vertical_deflection_end_2(3, parms);
parms.theta = ang_rotation(10, 3, parms);

total_ang_rotation(parms) 

function L = lift_force(alpha, parms)
% Calculate the lift force (N) on the whole tail wing, provided it is 
% oriented with an angle of attack of alpha (rad). 
L = 0.5*parms.rho_infty*parms.u_infty^2*parms.cl_alpha*alpha*parms.b*parms.c;
end 

function R = rotational_force(T, x)
% Consider the torque (N*m) exerted by the whole tail wing of length x (m) 
% to be experienced equally on both of its ends.  Calculate the rotational
% force due to the torque on either end of the beam.
R = (0.5*T)/(0.5*x);
end 

function FT1 = total_force_end_1(parms) 
% Calculate the total upwards force (N) exerted on the end of the beam, where
% the rotational force and the lift force add constructively (i.e., in the
% same direction). Suppose the total lift force exerted on the whole
% wing is experienced equally by both ends.
FT1 = parms.L/2 + parms.R; 
end 

function FT2 = total_force_end_2(parms)
% Calculate the total upwards force (N) exerted on the end of the beam, where
% the rotational force and the lift force add destructively (i.e., in 
% opposite directions). Suppose the total lift force exerted on the whole
% wing is experienced equally by both ends.
FT2 = parms.L/2 - parms.R; 
end 

function y1 = vertical_deflection_end_1(x, parms)
% Find the vertical deflection y (m) of end 1, given that it is subject to 
% a total vertical point load FT1 (N). The model is a cantilever beam of 
% length x, with one end fixed, the other end free and subject to a vertical 
% point load P. Assume the beam has an annular cross-section.
y1 = parms.FT1*x^3/(3*parms.E*parms.I); 
end 

function y2 = vertical_deflection_end_2(x, parms)
% Find the vertical deflection y (m) of end 2, given that it is subject to 
% a total vertical point load FT2 (N). The model is a cantilever beam of 
% length x, with one end fixed, the other end free and subject to a vertical 
% point load P. Assume the beam has an annular cross-section.
y2 = parms.FT2*x^3/(3*parms.E*parms.I); 
end 

function theta = ang_rotation(T, x, parms)
% Calculate the angular rotation (rad) of either end, as a result of the 
% torque T (N*m) exerted on the whole tail wing. 
theta = 0.5*T*x/(parms.G*parms.J);
end

function theta_tot = total_ang_rotation(parms)
% Calculate the total angular rotation (rad), experienced by the whole 
% tail wing.
theta_tot = parms.theta*2;
end 

function y_tot = total_vertical_deflection(parms)
y_tot = parms.y1 + parms.y2;
end 



