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
a_theta = A0(:,3);        % influencia de delta theta como "entrada"

%Para encontrar la funcion de transferencia
SI_A0 = s*eye(4) - A0;
Gp = simplify( C * inv(SI_A0)*a_theta );  % TF de la planta

% Separar numerador y denominador simbólicos (opcional)
% [numG, denG] = numden(Gp);
% numG = simplify(expand(numG));
% denG = simplify(expand(denG));

Gcb = Gp/(1+Gc*Gp);
Gcb = simplify(Gcb);
[~, denGcb] = numden(Gcb);
denGcb = simplify(expand(denGcb));

%% Diseño por amortiguamiento crítico
syms wn

% Polinomio dominante de segundo orden crítico
D_dom = s^2 + 2*wn*s + wn^2;

% Polos no dominantes (rápidos)
alpha = 5*wn;
beta = 6*wn;
D_des = expand(D_dom * (s + alpha) * (s + beta));

% Polinomio del lazo cerrado
D_lc = denGcb;

% % Igualación de coeficientes
% Pdiff = (expand(D_lc - D_des));
% coeffs_vec = coeffs(Pdiff, s);
% eqs = coeffs_vec == 0;

% Forma robusta: extraer polinomio en vector
[~,D_des_p] = numden(D_des);
poly_des = expand(D_des_p);
poly_lc  = expand(D_lc);

% Igualación por coeficientes (grado 4 -> 5 ecuaciones)
coeffs_lc = coeffs(poly_des,s,'All');
coeffs_des = coeffs(poly_des, s, 'All');

% alinear longitudes rellenando con ceros si es necesario
nL = length(coeffs_lc); nD = length(coeffs_des);
n = max(nL, nD);
coeffs_lc = [zeros(1,n-nL), coeffs_lc];
coeffs_des = [zeros(1,n-nD), coeffs_des];

% formar sistema de ecuaciones (vectorial)
eqs = coeffs_lc - coeffs_des;
eqs = simplify(eqs);           % ecuaciones simbólicas

% Resolver simbólicamente Kp, Kd, wn
sol = solve(eqs, [Kp Kd wn], 'Real', true, 'IgnoreAnalyticConstraints', true);

% recoger soluciones posibles
wn_sym = sol.wn;
Kp_sym = sol.Kp;
Kd_sym = sol.Kd;

% Kp = double(subs(Kp_sym, [wn alpha beta], [wn_val alpha_val beta_val]));
% Kd = double(subs(Kd_sym, [wn alpha beta], [wn_val alpha_val beta_val]));

% Si solve devuelve múltiples soluciones, filtrar:
found = false;
Kp = NaN; Kd = NaN; wn = NaN;
    for i = 1:length(wn_sym)
        try
            wn_i = double(wn_sym(i));
            Kp_i = double(Kp_sym(i));
            Kd_i = double(Kd_sym(i));
        catch
            continue;
        end
        if ~isreal(wn_i) || wn_i <= 0 || ~isreal(Kp_i) || ~isreal(Kd_i)
            continue;
        end
        % comprobar estabilidad del lazo cerrado
        poly_lc_sub = double(subs(poly_lc, {Kp,Kd,wn}, {Kp_i,Kd_i,wn_i}));
        % opcional: comprobar raices todas en LHP
        coeff_vec = sym2poly(poly_lc_sub);
        if isempty(coeff_vec)
            continue;
        end
        rts = roots(coeff_vec);
        if all(real(rts) < 0)
            Kp = Kp_i; Kd = Kd_i; wn = wn_i;
            found = true;
        end
    end

    if ~found
        warning('No se encontró solución válida para Kp,Kd,wn con los criterios dados.');
    end

end