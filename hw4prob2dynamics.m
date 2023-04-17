% Curran Robertson
% Lagrangian Equations of Motion of a Cart, Pendulum, Pulley System with
% No Dissipative Forces.
% April 16, 2023

clear all; clc; close all

% Symbols
syms m1 m2 m3 r theta dtheta x dx g l ddtheta ddx

% Constants

% Lagrangian ( L = T - V )
L = (1/2)*m1*r^2*dtheta^2 + (1/2)*m1*(dx+r*cos(theta))^2 + (1/2)*m1*(r*cos(theta))^2 + (1/2)*m2*dx^2 + (1/2)*m3*dx^2 - m1*g*r*cos(theta) - m3*g*(l-x);
L = expand(L);
L = simplify(L);
disp('L = ') 
pretty(L) 

% Partial Derivatives
    % DOF 1 : theta
d1 = diff(L, dtheta) ;
d2 = diff(L, theta);
    % DOF 2 : x
d4 = diff(L, dx)
d5 = diff(L, x)

% Time Derivatives
    % DOF 1 : theta
syms dtheta(t)
d3 = diff(subs(d1, dtheta, dtheta(t)), t);
    % DOF 2 : x
syms dx(t) theta(t)
d4_update = dx(t)*m1 + dx(t)*m2 + dx(t)*m3 + m1*r*cos(theta(t))
d6 = diff(d4_update, t)

% Eqns of Motion
    % DOF 1 : theta
eqn1 = d3 - d2 == 0;
    % DOF 2 : x
eqn2 = d6 - d5 == 0

% Solutions
    % DOF 1 : theta
ddtheta = solve(subs(eqn1, diff(dtheta(t), t), ddtheta), ddtheta);
disp('ddtheta = ')
pretty(ddtheta)
    % DOF 2 : x
ddx = solve(subs(eqn2, diff(dx(t), t), ddx), ddx);
disp('ddx = ')
pretty(ddx)
