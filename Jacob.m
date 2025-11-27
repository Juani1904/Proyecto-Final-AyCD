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

xt0 = 0;
w0 = 0;
Mt0 = M_x
ml0 = Mc_x;
bt0 = btw;
g0 = g;

A0 = subs(A, {xt,vt,th,w,Ftw,Mt,ml,l,bt,g}, {xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0});
B0 = subs(B, {xt,vt,th,w,Ftw,Mt,ml,l,bt,g}, {xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0});

% %% Buscamos la dinamica natural del sistema
% % Autovalores de la matriz A0
% lambda = eig(double(A0));
% 
% % Ordenar por real(lambda) (más cercano a 0 -> dominante)
% [~, idx] = sort(real(lambda), 'descend'); % mayor real -> más lento (menos negativo)
% lambda_ord = lambda(idx);
% 
% % Seleccionar el primer autovalor
% lambda_dom = lambda_ord(1);
% 
% % Si es complejo (par), buscar su conjugado y calcular wn_la
% if imag(lambda_dom) ~= 0
%     sigma = real(lambda_dom);
%     wd = abs(imag(lambda_dom));
%     wn_la = sqrt(sigma^2 + wd^2);
%     dseta = -sigma / wn_la;
% else
%     % polo real dominante
%     sigma = real(lambda_dom);   % negativo para estable
%     wn_la = abs(sigma);       % uso de |sigma| como medida de velocidad
%     dseta = 1;              % no aplica exactamente, placeholder
% end
% 
% %% Obtencion de Funcion de transferencia
% syms Kp Kd s real
% 
% Gc = Kp+Kd*s;
% 
% % Selección de salida y vector que representa la acción de delta-theta
% C = [0 1 0 0];            % salida y = delta v_t
% a_theta = A0(:,3);        % influencia de delta theta como "entrada"
% 
% %Para encontrar la funcion de transferencia
% SI_A0 = s*eye(4) - A0;
% Gp = simplify( C * inv(SI_A0)* B(:,1));  % TF de la planta
% 
% % Separar numerador y denominador simbólicos (opcional)
% % [numG, denG] = numden(Gp);
% % numG = simplify(expand(numG));
% % denG = simplify(expand(denG));
% 
% % 
% % Gct = (1.7024e+04)/(s^2)+(2.2754e+04)/s+(1.0849e+04);
% % Gct = simplify(Gct);
% % 
% % Gcb = Gp/(1+Gct*Gc*Gp); 
% 
% Gcb = Gp/(1+Gc*Gp)
% Gcb = simplify(Gcb);
% [~, denGcb] = numden(Gcb);
% denGcb = expand(denGcb);
% 
% %% Diseño por amortiguamiento crítico
% % Usamos una wn mayor a la de lazo abierto
% wn = 1.5*wn_la;
% 
% % Polinomio dominante de segundo orden crítico
% P_dom = s^2 + 2*wn*s + wn^2;
% 
% % Polos no dominantes (rápidos)
% alpha = 5*wn;
% beta = 6*wn;
% P_des = expand(P_dom * (s + alpha) * (s + beta));
% 
% % Polinomio del lazo cerrado
% P_lc = denGcb;
% 
% % % Igualación de coeficientes
% % Pdiff = (expand(D_lc - D_des));
% % coeffs_vec = coeffs(Pdiff, s);
% % eqs = coeffs_vec == 0;
% 
% % Forma robusta: extraer polinomio en vector
% % [~,D_des_p] = numden(D_des);
% % poly_des = expand(D_des);
% 
% % Igualación por coeficientes (grado 4 -> 5 ecuaciones)
% coeffs_lc = coeffs(P_lc,s,'All');
% coeffs_des = coeffs(P_des, s, 'All');
% 
% % alinear longitudes rellenando con ceros si es necesario
% nL = length(coeffs_lc); nD = length(coeffs_des);
% n = max(nL, nD);
% coeffs_lc = [zeros(1,n-nL), coeffs_lc];
% coeffs_des = [zeros(1,n-nD), coeffs_des];
% 
% % formar sistema de ecuaciones (vectorial)
% eqs = coeffs_lc - coeffs_des;
% eqs = simplify(eqs);           % ecuaciones simbólicas
% 
% % Resolver simbólicamente Kp, Kd
% sol = solve(eqs, [Kp Kd], 'Real', true);
% 
% % recoger soluciones posibles
% Kp_sol = sol.Kp;
% Kd_sol = sol.Kd;
% 
% % % Si solve devuelve múltiples soluciones, filtrar:
% % found = false;
% % Kp = NaN; Kd = NaN; wn = NaN;
% %     for i = 1:length(wn_sym)
% %         try
% %             wn_i = double(wn_sym(i));
% %             Kp_i = double(Kp_sym(i));
% %             Kd_i = double(Kd_sym(i));
% %         catch
% %             continue;
% %         end
% %         if ~isreal(wn_i) || wn_i <= 0 || ~isreal(Kp_i) || ~isreal(Kd_i)
% %             continue;
% %         end
% %         % comprobar estabilidad del lazo cerrado
% %         poly_lc_sub = double(subs(poly_lc, {Kp,Kd,wn}, {Kp_i,Kd_i,wn_i}));
% %         % opcional: comprobar raices todas en LHP
% %         coeff_vec = sym2poly(poly_lc_sub);
% %         if isempty(coeff_vec)
% %             continue;
% %         end
% %         rts = roots(coeff_vec);
% %         if all(real(rts) < 0)
% %             Kp = Kp_i; Kd = Kd_i; wn = wn_i;
% %             found = true;
% %         end
% %     end
% % 
% %     if ~found
% %         warning('No se encontró solución válida para Kp,Kd,wn con los criterios dados.');
% %     end
