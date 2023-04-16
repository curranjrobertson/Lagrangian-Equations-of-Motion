% Curran Robertson
% Lagrangian Equations of Motion of a Sliding Pendulum with Friction and
% Air Resistance
% April 16, 2023

clear all; clc; close all

% Symbols
syms m2 dfi fi x dx s ds t c mu ddfi dds ddx

% Constants
m1 = 1;
m2 = 1.5;
k = 25.6;
g = 9.81;
l = 0.8;

% Lagrangian ( L = T - V )
L = (1/2)*m1*dx^2 + (1/2)*m2*(dx+(l+s)*cos(fi)+ds*sin(fi))^2 + (1/2)*m2*((l+s)*sin(fi)+ds*cos(fi))^2 + (1/2)*m2*(l+s)^2*dfi^2 + m2*g*((l+s)*cos(fi)) - (1/2)*k*s^2;
L = expand(L);
L = simplify(L);
disp('L = ')
pretty(L)

% Partial Derivatives
    % DOF 1 : Fi
d1 = diff(L, dfi);
d2 = diff(L, fi);

    % DOF 2 : x
d4 = diff(L, dx);
d5 = diff(L, x);

    % DOF 3 : s
d7 = diff(L, ds);
d8 = diff(L, s);


% Time Derivatives
    % DOF 1 : Fi
syms dfi(t)
d3 = diff(subs(d1, dfi, dfi(t)), t);
    % DOF 2 : x
syms dx(t)
d6 = diff(subs(d4, dx, dx(t)), t);
    % DOF 3 : s
syms ds(t)
d9 = diff(subs(d7, ds, ds(t)), t);

% Eqns of Motion
    % DOF 1 : Fi
eqn1 = d3 - d2 == -c*dfi; % RHS = Dissipative force due to air resistance. Drag is proportional to velocity
    % DOF 2 : x
eqn2 = d6 - d5 == -mu*dx; % RHS = Dissipative force due to friction. Drag is proportional to velocity
    % DOF 3 : s
eqn3 = d9 - d8 == -c*ds; % RHS = Dissipative force due to air resistance. Drag is proportional to velocity

% Solutions
    % DOF 1 : Fi
ddfi = solve(subs(eqn1, diff(dfi(t), t), ddfi), ddfi);
disp('ddfi = ')
pretty(ddfi)
    % DOF 2 : x
ddx = solve(subs(eqn2, diff(dx(t), t), ddx), ddx);
disp('ddx = ')
pretty(ddx)
    % DOF 3 : s
dds = solve(subs(eqn3, diff(ds(t), t), dds), dds);
disp('dds = ')
pretty(dds)
