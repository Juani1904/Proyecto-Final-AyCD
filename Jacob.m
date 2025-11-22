function [Kp, Kd, wn] = jacob(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)

syms xt vt th w Ftw real
syms Mt ml l bt g real

%%Linealizacion Jacobiana 

% Estados y entrada
x = [xt; vt; th; w];
u = Ftw;

% Definiciones auxiliares
Dx  = Mt + ml*sin(th)^2;
Nx = Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th);
Nth = -(Mt+ml*sin(th)^2)*g*sin(th)-cos(th)*(Ftw-bt*vt+ml*l*w^2*sin(th)+ml*g*sin(th)*cos(th));
Dth = l*(Mt+ml*sin(th)^2);

% Dinámica
f1 = vt;
f2 = Nx / Dx;
f3 = w;
f4 = Nth / Dth;

f = [f1; f2; f3; f4];

% Jacobianos
A = jacobian(f, x);
B = jacobian(f, u);

% (opcional) simplificación
A = simplify(A);
B = simplify(B);

% Modelo LPV evaluado en el punto de operacion
A0 = subs(A, {xt,vt,th,w,Ftw,Mt,ml,l,bt,g}, {xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0});
B0 = subs(B, {xt,vt,th,w,Ftw,Mt,ml,l,bt,g}, {xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0});


%% Obtencion de Funcion de transferencia
syms Kp Kd s real

Gc = Kp+Kd*s;

% Selección de salida y vector que representa la acción de delta-theta
C = [0 1 0 0];            % salida y = delta v_t
a_theta = A(:,3);        % influencia de delta theta como "entrada"

%Para encontrar la funcion de transferencia
SI_A0 = s*eye(4) - A0;
Gp = simplify( C * inv(SI_A0) * a_theta );  % TF de la planta

% Separar numerador y denominador simbólicos (opcional)
[numG, denG] = numden(Gp);
numG = simplify(expand(numG));
denG = simplify(expand(denG));

Gcb = Gp/(1+Gc*Gp);
Gcb = simplify(Gcb);
[numGcb, denGcb] = numden(Gcb);
numGcb = simplify(expand(numGcb));
denGcb = simplify(expand(denGcb));

%% Diseño por amortiguamiento crítico
syms wn alpha beta

% Polinomio dominante de segundo orden crítico
D_dom = s^2 + 2*wn*s + wn^2;

% Polos no dominantes (rápidos)
D_des = expand(D_dom * (s + alpha) * (s + beta));

% Polinomio del lazo cerrado
D_lc = denGcb;

% Igualación de coeficientes
eqs = coeffs(expand(D_lc - D_des), s) == 0;

% Resolver simbólicamente Kp y Kd
sol = solve(eqs, [Kp Kd], 'Real', true);

Kp_sym = sol.Kp;
Kd_sym = sol.Kd;

%% Evaluación numérica
% Definición de parámetros deseados
wn_val    = 4;               % ejemplo
alpha_val = 5*wn_val;
beta_val  = 6*wn_val;

Kp = double(subs(Kp_sym, [wn alpha beta], [wn_val alpha_val beta_val]));
Kd = double(subs(Kd_sym, [wn alpha beta], [wn_val alpha_val beta_val]));

wn = wn_val;   % se devuelve para referencia

end