function [A, B] = LIN_JACOB(Mt, ml, l, bt, g, xdot, xddot)
% LIN_JACOB
%   Calcula las matrices A y B linealizadas para un sistema carro-péndulo colgante,
%   considerando un punto de equilibrio con ángulo theta_eq y velocidad del carro xdot_eq.
%
%   Parámetros de entrada:
%       Mt        -> masa del carro [kg]
%       ml        -> masa del péndulo [kg]
%       l         -> longitud del péndulo [m]
%       bt        -> coeficiente de fricción viscosa del carro [N·s/m]
%       g         -> gravedad [m/s^2]
%       Ft        -> Fuerza de tensión del cable
%
%   Salidas:
%       A, B  -> matrices linealizadas en el punto de equilibrio
%
%   Ecuaciones base:
%       (Mt + ml) * x_ddot + ml*l*(theta_ddot*cos(theta) - theta_dot^2*sin(theta)) = Ft - bt*x_dot
%       ml*l^2*theta_ddot + ml*l*x_ddot*cos(theta) + ml*g*l*sin(theta) = 0
%
%

%% === Punto de equilibrio dinámico ===
%xddot = -g * tan(theta_eq);           % aceleración necesaria para mantener ángulo constante
theta_eq = atan2(xddot/-g);             %Calculo de theta de equilibrio para cierta aceleracion

%% === Seno y coseno ===
c = cos(theta_eq);
s = sin(theta_eq);

%% === Cálculo de x'' ===
% === Calculo del Numerador ===
Nx = Ft - bt*xdot_eq + m*g*s*c;
dNx = m*g*(c^2 - s^2);
% === Cálculo de denominador común ===
Dx = Mt + m*s^2;
dDx = m*2*s*c;

%% === Cálculo de theta'' ===
% === Cálculo del Numerador ===
Nth = -(Mt+m)*g*s - c*(Ft - bt*v);
dNth = -(Mt+m)*g*c + s*(Ft - bt*v);
% === Cálculo del Denominador ===
Dth = l*(Mt + m*s^2);
dDth = l*(2*m*s*c);

%% === Coeficientes de matrices ===
a22 = -bt / Dx;
a23 = (dNx*Dx - Nx*dDx) / Dx^2;
b2 = 1 / (Mt + m*s^2);

a42 = (c*bt) / Dth;
a43 = (dNth*Dth - Nth*dDth) / (Dth^2);
b4 = -c / (l*(Mt + m*s^2));


%% === Matrices A y B ===
A = [  0      1      0      0;
       0    a22      a23    0;
       0      0      0      1;
       0    a42      a43    0];

B = [ 0;
      b2;
      0;
      b4];

end