function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym, Gcb_sym] = ...
    jacob_gs(xt0,vt0,at0,w0,Mt0,ml0,l0,bt0,g0)

th0 = atan2(-at0,g0);

Ftw0 = (Mt0 + ml0)*at0 + bt0*vt0;

syms xt vt th w Ftw Mt ml l bt g dd_xt dd_th s real

%% Estados e input
x = [xt; vt; th; w]; %Estados x1=xt x2=d_xt=vt x3=th x4=d_th=w
u = Ftw;

%% Modelo no lineal exacto - Modelo dinámico
eq1 = (Mt + ml)*dd_xt + ml*l*dd_th*cos(th) - ml*l*w^2*sin(th) == Ftw - bt*vt;
eq2 = l*dd_th + dd_xt*cos(th) + g*sin(th) == 0;

sol = solve([eq1, eq2], [dd_xt, dd_th]);

dd_xt_expr  = simplify(sol.dd_xt);
dd_th_expr = simplify(sol.dd_th);

%% Dinamica de estados: xdot = f(x,u) MODELO EN ESPACIO DE ESTADOS
f = [vt; dd_xt_expr; w; dd_th_expr];

%% Jacobianos
A = simplify(jacobian(f, x));
B = simplify(jacobian(f, u));

%% Punto de operacion
%solF = solve(eq2,Ftw);
%Ftw_expr = simplify(solF.Ftw);
%Ftw0= subs(Ftw_expr,[])

vars_sym = [xt; vt; th; w; Ftw; Mt; ml; l; bt; g];
vals_num = [xt0; vt0; th0; w0; Ftw0; Mt0; ml0; l0; bt0; g0];

A0 = subs(A, vars_sym, vals_num);
B0 = subs(B, vars_sym, vals_num);

A0_num = double(A0);
B0_num = double(B0);

%% Planta local linealizada Gp(s) = Theta(s)/Vt(s)
Cth = [0 0 1 0];
Cv  = [0 1 0 0];

GthF = minreal(tf(ss(A0_num, B0_num, Cth, 0)));
GvF  = minreal(tf(ss(A0_num, B0_num, Cv, 0)));

Gp_tf = minreal(GthF / GvF);

[numGp, denGp] = tfdata(Gp_tf, 'v');
numGp = quitar_ceros_inicio(real(numGp));
denGp = quitar_ceros_inicio(real(denGp));

if numel(denGp) ~= 3 || numel(numGp) ~= 2
    error('La sintonia serie implementada requiere Gp(s) de orden 2 con numerador de orden 1.');
end

Gp_sym = poly2sym(numGp, s) / poly2sym(denGp, s);

%% Coeficientes de Gp(s) = (b1*s + b0)/(a2*s^2 + a1*s + a0)
a2 = denGp(1);
a1 = denGp(2);
a0 = denGp(3);

b1 = numGp(1);
b0 = numGp(2);

%% Parametros de sintonia serie sobre Gcb(s)
% Gcb(s) = Gp(s) / (1 + Gc(s)*Gp(s))
% Gc(s)  = Kd*s + Kp
%
% Denominador de Gcb:
% Dcb(s) = Dp(s) + Gc(s)Np(s)
%        = (a2 + Kd*b1)s^2 + (a1 + Kd*b0 + Kp*b1)s + (a0 + Kp*b0)
%
% Se impone el polinomio deseado normalizado:
% Pdes(s) = s^2 + alpha*s + beta
% con alpha = n_p*w_pos_p y beta = w_pos_p^2

polos_gp = roots(denGp);
polos_gp = polos_gp(isfinite(polos_gp));

if isempty(polos_gp)
    w_nat_gp = sqrt(abs(a0/a2));
else
    w_nat_gp = max(abs(polos_gp));
end

lambda_ss = 10;
w_pos_p   = lambda_ss * w_nat_gp;
n_p       = 1.92;

alpha = n_p * w_pos_p;
beta  = w_pos_p^2;

%% Calculo de Kp y Kd igualando el denominador de Gcb al polinomio deseado
% Ecuaciones:
% (a1 + Kd*b0 + Kp*b1)/(a2 + Kd*b1) = alpha
% (a0 + Kp*b0)/(a2 + Kd*b1)         = beta

M = [b1,            (b0 - alpha*b1);
     b0,            (-beta*b1)];

rhs = [alpha*a2 - a1;
       beta*a2  - a0];

if abs(det(M)) < 1e-12 || rcond(M) < 1e-12
    error('El sistema para calcular Kp y Kd esta mal condicionado.');
end

solK = M \ rhs;

Kp_num = real(solK(1));
Kd_num = real(solK(2));

if ~isfinite(Kp_num) || ~isfinite(Kd_num)
    error('No se pudieron obtener ganancias finitas.');
end

%% Controlador PD
Gc_num = [Kd_num Kp_num];

%% Funcion de transferencia a lazo cerrado
num_cb = numGp;
den_cb = suma_poly(denGp, conv(Gc_num, numGp));

if isempty(den_cb) || abs(den_cb(1)) < 1e-12
    error('El denominador de lazo cerrado es invalido.');
end

den_cb = real(den_cb / den_cb(1));
num_cb = real(num_cb / den_cb(1)); %#ok<NASGU>

Plc_sym  = poly2sym(den_cb, s);
Gcb_sym  = poly2sym(numGp, s) / poly2sym(den_cb, s);

%% Polinomio deseado
Pdes = [1 alpha beta];
Pdes_sym = poly2sym(Pdes, s);

%% Chequeo de estabilidad
r_cl = roots(den_cb);
if any(real(r_cl) >= 0)
    warning('El lazo cerrado calculado no quedo estable para este punto de operacion.');
end

end

function p = quitar_ceros_inicio(p)
tol = 1e-12;
idx = find(abs(p) > tol, 1, 'first');
if isempty(idx)
    p = 0;
else
    p = p(idx:end);
end
end

function p = suma_poly(p1, p2)
n = max(length(p1), length(p2));
p1 = [zeros(1, n-length(p1)), p1];
p2 = [zeros(1, n-length(p2)), p2];
p  = quitar_ceros_inicio(real(p1 + p2));
end


% 
% function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym] = ...
%     jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)
% 
% syms xt vt th w Ftw real
% syms Mt ml l bt g real
% syms xdd thdd real
% syms s Kp Kd real
% 
% %% Estados e input
% x = [xt; vt; th; w];
% u = Ftw;
% 
% %% Modelo no lineal exacto
% eq1 = (Mt + ml)*xdd + ml*l*thdd*cos(th) - ml*l*w^2*sin(th) == Ftw - bt*vt;
% eq2 = l*thdd + xdd*cos(th) + g*sin(th) == 0;
% 
% sol = solve([eq1, eq2], [xdd, thdd]);
% 
% xdd_expr  = simplify(sol.xdd);
% thdd_expr = simplify(sol.thdd);
% 
% disp(xdd_expr);
% disp(thdd_expr);
% 
% f1 = vt;
% f2 = xdd_expr;
% f3 = w;
% f4 = thdd_expr;
% 
% f = [f1; f2; f3; f4];
% 
% %% Jacobianos
% A = simplify(jacobian(f, x));
% B = simplify(jacobian(f, u));
% 
% %% Punto de operación
% vars_sym = {xt, vt, th, w, Ftw, Mt, ml, l, bt, g};
% vals_num = {xt0, vt0, th0, w0, Ftw0, Mt0, ml0, l0, bt0, g0};
% 
% A0 = simplify(subs(A, vars_sym, vals_num));
% B0 = simplify(subs(B, vars_sym, vals_num));
% 
% A0_num = double(A0);
% B0_num = double(B0);
% 
% %% Planta local linealizada
% sys_lin = ss(A0_num, B0_num, eye(4), zeros(4,1));
% 
% % Theta/Ftw
% Cth = [0 0 1 0];
% GthF = minreal(tf(ss(A0_num, B0_num, Cth, 0)));
% 
% % Vt/Ftw
% Cv  = [0 1 0 0];
% GvF = minreal(tf(ss(A0_num, B0_num, Cv, 0)));
% 
% % Planta efectiva para el anti-sway: Theta/Vt
% Gp_tf = minreal(GthF / GvF);
% 
% [numGp, denGp] = tfdata(Gp_tf, 'v');
% numGp = quitar_ceros_inicio(real(numGp));
% denGp = quitar_ceros_inicio(real(denGp));
% 
% Gp_sym = poly2sym(numGp, s) / poly2sym(denGp, s);
% 
% %% Polinomio deseado
% n = length(denGp) - 1;
% 
% if n < 2
%     warning('La planta local Theta/Vt tiene orden menor que 2.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Pdes_sym = sym(0);
%     Plc_sym  = sym(0);
%     return
% end
% 
% p_ol = pole(Gp_tf);
% rp = abs(real(p_ol));
% rp = rp(rp > 1e-4);
% 
% if isempty(rp)
%     wn_base = 0.5;
% else
%     wn_base = min(rp);
% end
% 
% wn   = max(1.5*wn_base, 0.5);
% zeta = 0.9;
% 
% sigma = zeta*wn;
% wd    = wn*sqrt(max(1-zeta^2,0));
% 
% if wd > 1e-8
%     p_dom = [-sigma + 1i*wd, -sigma - 1i*wd];
% else
%     p_dom = [-wn, -1.2*wn];
% end
% 
% p_extra = zeros(1, max(n-2,0));
% for k = 1:length(p_extra)
%     p_extra(k) = -(3+k)*wn;
% end
% 
% p_des = [p_dom, p_extra];
% Pdes  = real(poly(p_des));
% Pdes  = Pdes / Pdes(1);
% Pdes_sym = poly2sym(Pdes, s);
% 
% %% Búsqueda numérica de Kp y Kd sobre el lazo cerrado
% starts = [...
%      0.1   0.01;
%      1     0.1;
%      10    1;
%      100   10;
%     -0.1   0.01;
%     -1     0.1;
%     -10    1;
%      0.1  -0.01;
%      1    -0.1;
%      10   -1;
%     -0.1  -0.01;
%     -1    -0.1;
%     -10   -1];
% 
% bestJ   = inf;
% bestKp  = NaN;
% bestKd  = NaN;
% bestDen = [];
% 
% opts = optimset('Display','off', ...
%                 'TolX',1e-8, ...
%                 'TolFun',1e-8, ...
%                 'MaxIter',500, ...
%                 'MaxFunEvals',3000);
% 
% for ii = 1:size(starts,1)
%     x0 = starts(ii,:);
% 
%     try
%         xopt = fminsearch(@(x)costo_pd_cl(x, numGp, denGp, Pdes), x0, opts);
% 
%         Kp_i = xopt(1);
%         Kd_i = xopt(2);
% 
%         den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp_i, Kd_i);
%         rts = roots(den_cl);
% 
%         if all(real(rts) < -1e-5)
%             J = costo_pd_cl([Kp_i Kd_i], numGp, denGp, Pdes);
%             if J < bestJ
%                 bestJ   = J;
%                 bestKp  = Kp_i;
%                 bestKd  = Kd_i;
%                 bestDen = den_cl;
%             end
%         end
%     catch
%     end
% end
% 
% if isfinite(bestJ)
%     Kp_num = bestKp;
%     Kd_num = bestKd;
%     Plc_sym = poly2sym(bestDen, s);
% else
%     warning('No se encontró una solución estable para Kp y Kd.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Plc_sym = sym(0);
% end
% 
% end
% 
% function J = costo_pd_cl(x, numGp, denGp, Pdes)
% Kp = x(1);
% Kd = x(2);
% 
% if any(~isfinite([Kp Kd]))
%     J = 1e20;
%     return
% end
% 
% den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp, Kd);
% 
% if isempty(den_cl) || any(~isfinite(den_cl)) || abs(den_cl(1)) < 1e-12
%     J = 1e20;
%     return
% end
% 
% den_cl = real(den_cl / den_cl(1));
% 
% n = max(length(den_cl), length(Pdes));
% den_cl_pad = [zeros(1,n-length(den_cl)), den_cl];
% Pdes_pad   = [zeros(1,n-length(Pdes)),   Pdes];
% 
% err_coeff = norm(den_cl_pad - Pdes_pad, 2)^2;
% 
% rts = roots(den_cl);
% pen_est  = 1e6 * sum(max(real(rts)+1e-5,0).^2);
% pen_kd   = 1e-3 / (abs(Kd) + 1e-4);
% pen_gain = 1e-8 * (Kp^2 + Kd^2);
% 
% J = err_coeff + pen_est + pen_kd + pen_gain;
% end
% 
% function den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp, Kd)
% Gc_num = [Kd Kp];
% 
% num_ol = conv(Gc_num, numGp);
% den_ol = denGp;
% 
% n = max(length(num_ol), length(den_ol));
% num_ol = [zeros(1,n-length(num_ol)), num_ol];
% den_ol = [zeros(1,n-length(den_ol)), den_ol];
% 
% den_cl = den_ol + num_ol;
% den_cl = quitar_ceros_inicio(real(den_cl));
% end
% 
% function p = quitar_ceros_inicio(p)
% idx = find(abs(p) > 1e-12, 1, 'first');
% if isempty(idx)
%     p = 0;
% else
%     p = p(idx:end);
% end
% end
% function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym] = ...
%     jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)
% 
% syms xt vt th w Ftw real
% syms Mt ml l bt g real
% syms xdd thdd real
% syms s real
% 
% %% Estados e input
% x = [xt; vt; th; w];
% u = Ftw;
% 
% %% Ecuaciones no lineales exactas del sistema
% eq1 = (Mt + ml)*xdd + ml*l*thdd*cos(th) - ml*l*w^2*sin(th) == Ftw - bt*vt;
% eq2 = l*thdd + xdd*cos(th) + g*sin(th) == 0;
% 
% %% Resolver simbólicamente para xdd y thdd
% sol = solve([eq1, eq2], [xdd, thdd]);
% 
% xdd_expr  = simplify(sol.xdd);
% thdd_expr = simplify(sol.thdd);
% 
% %% Modelo no lineal en espacio de estados
% f1 = vt;
% f2 = xdd_expr;
% f3 = w;
% f4 = thdd_expr;
% 
% f = [f1; f2; f3; f4];
% 
% %% Jacobianos exactos
% A = simplify(jacobian(f, x));
% B = simplify(jacobian(f, u));
% 
% %% Evaluación en punto de operación
% vars_sym = {xt, vt, th, w, Ftw, Mt, ml, l, bt, g};
% vals_num = {xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0};
% 
% A0 = simplify(subs(A, vars_sym, vals_num));
% B0 = simplify(subs(B, vars_sym, vals_num));
% 
% A0_num = double(A0);
% B0_num = double(B0);
% 
% %% Salida angular
% C0 = [0 1 0 0];
% D0 = 0;
% 
% %% Planta local Theta(s)/Ftw(s)
% Gp_ss = ss(A0_num, B0_num, C0, D0);
% Gp_tf = minreal(tf(Gp_ss));
% 
% [numGp, denGp] = tfdata(Gp_tf, 'v');
% numGp = quitar_ceros_inicio(real(numGp));
% denGp = quitar_ceros_inicio(real(denGp));
% 
% Gp_sym = poly2sym(numGp, s) / poly2sym(denGp, s);
% 
% %% Polinomio deseado
% n = length(denGp) - 1;
% 
% if n < 2
%     warning('La planta local tiene orden menor que 2.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Pdes_sym = sym(0);
%     Plc_sym = sym(0);
%     return
% end
% 
% lambda = eig(A0_num);
% real_part = abs(real(lambda));
% real_part = real_part(real_part > 1e-4);
% 
% if isempty(real_part)
%     wn_base = 0.5;
% else
%     wn_base = min(real_part);
% end
% 
% wn   = max(1.5*wn_base, 0.5);
% zeta = 0.9;
% 
% sigma = zeta*wn;
% wd    = wn*sqrt(max(1-zeta^2,0));
% 
% if wd > 1e-8
%     p_dom = [-sigma + 1i*wd, -sigma - 1i*wd];
% else
%     p_dom = [-wn, -1.2*wn];
% end
% 
% p_extra = zeros(1, max(n-2,0));
% for k = 1:length(p_extra)
%     p_extra(k) = -(3+k)*wn;
% end
% 
% p_des = [p_dom, p_extra];
% Pdes  = real(poly(p_des));
% Pdes  = Pdes / Pdes(1);
% Pdes_sym = poly2sym(Pdes, s);
% 
% %% Búsqueda numérica de Kp y Kd para Gc(s)=Kp+Kd*s
% starts = [...
%      0.1   0.01;
%      1     0.1;
%      10    1;
%      100   10;
%     -0.1   0.01;
%     -1     0.1;
%     -10    1;
%      0.1  -0.01;
%      1    -0.1;
%      10   -1;
%     -0.1  -0.01;
%     -1    -0.1;
%     -10   -1];
% 
% bestJ   = inf;
% bestKp  = NaN;
% bestKd  = NaN;
% bestDen = [];
% 
% opts = optimset('Display','off', ...
%                 'TolX',1e-8, ...
%                 'TolFun',1e-8, ...
%                 'MaxIter',500, ...
%                 'MaxFunEvals',3000);
% 
% for ii = 1:size(starts,1)
%     x0 = starts(ii,:);
% 
%     try
%         xopt = fminsearch(@(x)costo_pd(x, numGp, denGp, Pdes), x0, opts);
% 
%         Kp_i = xopt(1);
%         Kd_i = xopt(2);
% 
%         den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp_i, Kd_i);
%         rts = roots(den_cl);
% 
%         if all(real(rts) < -1e-5)
%             J = costo_pd([Kp_i Kd_i], numGp, denGp, Pdes);
% 
%             if J < bestJ
%                 bestJ   = J;
%                 bestKp  = Kp_i;
%                 bestKd  = Kd_i;
%                 bestDen = den_cl;
%             end
%         end
%     catch
%     end
% end
% 
% if isfinite(bestJ)
%     Kp_num = bestKp;
%     Kd_num = bestKd;
%     Plc_sym = poly2sym(bestDen, s);
% else
%     warning('No se encontró una solución estable para Kp y Kd.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Plc_sym = sym(0);
% end
% 
% end
% 
% function J = costo_pd(x, numGp, denGp, Pdes)
% Kp = x(1);
% Kd = x(2);
% 
% if any(~isfinite([Kp Kd]))
%     J = 1e20;
%     return
% end
% 
% den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp, Kd);
% 
% if isempty(den_cl) || any(~isfinite(den_cl)) || abs(den_cl(1)) < 1e-12
%     J = 1e20;
%     return
% end
% 
% den_cl = real(den_cl / den_cl(1));
% 
% n = max(length(den_cl), length(Pdes));
% den_cl_pad = [zeros(1,n-length(den_cl)), den_cl];
% Pdes_pad   = [zeros(1,n-length(Pdes)),   Pdes];
% 
% err_coeff = norm(den_cl_pad - Pdes_pad, 2)^2;
% 
% rts = roots(den_cl);
% pen_est  = 1e6 * sum(max(real(rts)+1e-5,0).^2);
% pen_kd   = 1e-3 / (abs(Kd) + 1e-4);
% pen_gain = 1e-8 * (Kp^2 + Kd^2);
% 
% J = err_coeff + pen_est + pen_kd + pen_gain;
% end
% 
% function den_cl = den_lazo_cerrado_pd(numGp, denGp, Kp, Kd)
% Gc_num = [Kd Kp];
% 
% num_ol = conv(Gc_num, numGp);
% den_ol = denGp;
% 
% n = max(length(num_ol), length(den_ol));
% num_ol = [zeros(1,n-length(num_ol)), num_ol];
% den_ol = [zeros(1,n-length(den_ol)), den_ol];
% 
% den_cl = den_ol + num_ol;
% den_cl = quitar_ceros_inicio(real(den_cl));
% end
% 
% function p = quitar_ceros_inicio(p)
% idx = find(abs(p) > 1e-12, 1, 'first');
% if isempty(idx)
%     p = 0;
% else
%     p = p(idx:end);
% end
% end

% function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym] = ...
%     jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)
% % Ftw0 se mantiene solo por compatibilidad con el pipeline existente
% 
% % Se asume que estas variables existen en base workspace:
% % it, rtd, btm, btd, Jtm_tb, Jtd
% % Opcional: Ttb0
% 
% it0     = evalin('base','it');
% rtd0    = evalin('base','rtd');
% btm0    = evalin('base','btm');
% btd0    = evalin('base','btd');
% Jtm_tb0 = evalin('base','Jtm_tb');
% Jtd0    = evalin('base','Jtd');
% 
% if evalin('base','exist(''Ttb0'',''var'')')
%     Ttb0 = evalin('base','Ttb0');
% else
%     Ttb0 = 0;
% end
% 
% syms xt vt th w Ttm real
% syms Mt ml l bt g real
% syms it rtd btm btd Jtm_tb Jtd Ttb real
% syms s real
% 
% %% Estados e input real
% x = [xt; vt; th; w];
% u = Ttm;
% 
% %% Parámetros equivalentes del accionamiento
% Aeq = it/rtd;
% Beq = (it^2*btm + btd)/rtd^2;
% Meq = (it^2*Jtm_tb + Jtd)/rtd^2;
% 
% %% Modelo no lineal extendido con entrada Ttm
% Den = Mt + ml*sin(th)^2 + Meq;
% 
% Phi = ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th);
% 
% f1 = vt;
% f2 = ( Aeq*(Ttm + Ttb) - (bt + Beq)*vt + Phi ) / Den;
% f3 = w;
% f4 = -( g*sin(th) + cos(th)*f2 ) / l;
% 
% f = [f1; f2; f3; f4];
% 
% %% Jacobianos
% A = simplify(jacobian(f,x));
% B = simplify(jacobian(f,u));
% 
% %% Punto de operación
% vars_sym = {xt,vt,th,w,Ttm,Mt,ml,l,bt,g,it,rtd,btm,btd,Jtm_tb,Jtd,Ttb};
% vals_num = {xt0,vt0,th0,w0,0,Mt0,ml0,l0,bt0,g0,it0,rtd0,btm0,btd0,Jtm_tb0,Jtd0,Ttb0};
% 
% A0 = simplify(subs(A, vars_sym, vals_num));
% B0 = simplify(subs(B, vars_sym, vals_num));
% 
% A0_num = double(A0);
% B0_num = double(B0);
% 
% %% Salida angular
% C0 = [0 0 1 0];
% D0 = 0;
% 
% %% Planta local Theta(s)/Ttm(s)
% Gp_ss = ss(A0_num,B0_num,C0,D0);
% Gp_tf = minreal(tf(Gp_ss));
% 
% [numGp, denGp] = tfdata(Gp_tf,'v');
% numGp = quitar_ceros_inicio(real(numGp));
% denGp = quitar_ceros_inicio(real(denGp));
% 
% Gp_sym = poly2sym(numGp,s) / poly2sym(denGp,s);
% 
% %% Polinomio deseado
% n = length(denGp) - 1;
% 
% if n < 2
%     warning('La planta local tiene orden menor que 2.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Pdes_sym = sym(0);
%     Plc_sym  = sym(0);
%     return
% end
% 
% lambda = eig(A0_num);
% real_part = abs(real(lambda));
% real_part = real_part(real_part > 1e-4);
% 
% if isempty(real_part)
%     wn_base = 0.5;
% else
%     wn_base = min(real_part);
% end
% 
% wn   = max(1.5*wn_base, 0.5);
% zeta = 0.9;
% 
% sigma = zeta*wn;
% wd    = wn*sqrt(max(1-zeta^2,0));
% 
% if wd > 1e-8
%     p_dom = [-sigma + 1i*wd, -sigma - 1i*wd];
% else
%     p_dom = [-wn, -1.2*wn];
% end
% 
% p_extra = zeros(1,max(n-2,0));
% for k = 1:length(p_extra)
%     p_extra(k) = -(3+k)*wn;
% end
% 
% p_des = [p_dom, p_extra];
% Pdes  = real(poly(p_des));
% Pdes  = Pdes / Pdes(1);
% Pdes_sym = poly2sym(Pdes,s);
% 
% %% Búsqueda numérica de Kp y Kd
% starts = [...
%      0.1   0.01;
%      1     0.1;
%      10    1;
%      100   10;
%     -0.1   0.01;
%     -1     0.1;
%     -10    1;
%      0.1  -0.01;
%      1    -0.1;
%      10   -1;
%     -0.1  -0.01;
%     -1    -0.1;
%     -10   -1];
% 
% bestJ   = inf;
% bestKp  = NaN;
% bestKd  = NaN;
% bestDen = [];
% 
% opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',500,'MaxFunEvals',3000);
% 
% for ii = 1:size(starts,1)
%     x0 = starts(ii,:);
% 
%     try
%         xopt = fminsearch(@(x)costo_pd(x,numGp,denGp,Pdes),x0,opts);
%         Kp_i = xopt(1);
%         Kd_i = xopt(2);
% 
%         den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp_i,Kd_i);
%         rts = roots(den_cl);
% 
%         if all(real(rts) < -1e-5)
%             J = costo_pd([Kp_i Kd_i],numGp,denGp,Pdes);
%             if J < bestJ
%                 bestJ   = J;
%                 bestKp  = Kp_i;
%                 bestKd  = Kd_i;
%                 bestDen = den_cl;
%             end
%         end
%     catch
%     end
% end
% 
% if isfinite(bestJ)
%     Kp_num = bestKp;
%     Kd_num = bestKd;
%     Plc_sym = poly2sym(bestDen,s);
% else
%     warning('No se encontró una solución estable para Kp y Kd.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Plc_sym = sym(0);
% end
% 
% end
% 
% function J = costo_pd(x,numGp,denGp,Pdes)
% Kp = x(1);
% Kd = x(2);
% 
% if any(~isfinite([Kp Kd]))
%     J = 1e20;
%     return
% end
% 
% den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd);
% 
% if isempty(den_cl) || any(~isfinite(den_cl)) || abs(den_cl(1)) < 1e-12
%     J = 1e20;
%     return
% end
% 
% den_cl = real(den_cl / den_cl(1));
% 
% n = max(length(den_cl),length(Pdes));
% den_cl_pad = [zeros(1,n-length(den_cl)), den_cl];
% Pdes_pad   = [zeros(1,n-length(Pdes)),   Pdes];
% 
% err_coeff = norm(den_cl_pad - Pdes_pad,2)^2;
% 
% rts = roots(den_cl);
% pen_est  = 1e6 * sum(max(real(rts)+1e-5,0).^2);
% pen_kd   = 1e-3 / (abs(Kd) + 1e-4);
% pen_gain = 1e-8 * (Kp^2 + Kd^2);
% 
% J = err_coeff + pen_est + pen_kd + pen_gain;
% end
% 
% function den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd)
% Gc_num = [Kd Kp];
% 
% num_ol = conv(Gc_num,numGp);
% den_ol = denGp;
% 
% n = max(length(num_ol),length(den_ol));
% num_ol = [zeros(1,n-length(num_ol)), num_ol];
% den_ol = [zeros(1,n-length(den_ol)), den_ol];
% 
% den_cl = den_ol + num_ol;
% den_cl = quitar_ceros_inicio(real(den_cl));
% end
% 
% function p = quitar_ceros_inicio(p)
% idx = find(abs(p) > 1e-12,1,'first');
% if isempty(idx)
%     p = 0;
% else
%     p = p(idx:end);
% end
% end



% function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym] = ...
%     jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)
% % Ftw0 se mantiene solo por compatibilidad
% 
% syms xt vt th w Tm_a Tm real
% syms Mt ml l bt g real
% syms s real
% syms tau_a K_a R_eq eta_eq real
% 
% x = [xt; vt; th; w; Tm_a];
% u = Tm;
% 
% Tm0   = 0;
% Tm_a0 = 0;
% 
% f_a = (-Tm_a + K_a*Tm)/tau_a;
% Ftw = eta_eq*Tm_a/R_eq;
% 
% Dx  = Mt + ml*sin(th)^2;
% Nx  = Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th);
% 
% Nth = -(Mt + ml*sin(th)^2)*g*sin(th) ...
%       - cos(th)*(Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th));
% 
% Dth = l*(Mt + ml*sin(th)^2);
% 
% f1 = vt;
% f2 = Nx / Dx;
% f3 = w;
% f4 = Nth / Dth;
% 
% f = [f1; f2; f3; f4; f_a];
% 
% A = simplify(jacobian(f,x));
% B = simplify(jacobian(f,u));
% 
% vars_sym = {xt,vt,th,w,Tm_a,Tm,Mt,ml,l,bt,g,tau_a,K_a,R_eq,eta_eq};
% vals_num = {xt0,vt0,th0,w0,Tm_a0,Tm0,Mt0,ml0,l0,bt0,g0,0.05,1,0.5,1};
% 
% A0 = simplify(subs(A, vars_sym, vals_num));
% B0 = simplify(subs(B, vars_sym, vals_num));
% 
% A0_num = double(A0);
% B0_num = double(B0);
% 
% C0 = [0 0 1 0 0];
% D0 = 0;
% 
% Gp_ss = ss(A0_num,B0_num,C0,D0);
% Gp_tf = minreal(tf(Gp_ss));
% 
% [numGp, denGp] = tfdata(Gp_tf,'v');
% numGp = quitar_ceros_inicio(real(numGp));
% denGp = quitar_ceros_inicio(real(denGp));
% 
% Gp_sym = poly2sym(numGp,s) / poly2sym(denGp,s);
% 
% n = length(denGp) - 1;
% 
% if n < 2
%     warning('La planta local tiene orden menor que 2.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Pdes_sym = sym(0);
%     Plc_sym  = sym(0);
%     return
% end
% 
% lambda = eig(A0_num);
% real_part = abs(real(lambda));
% real_part = real_part(real_part > 1e-4);
% 
% if isempty(real_part)
%     wn_base = 0.5;
% else
%     wn_base = min(real_part);
% end
% 
% wn   = max(1.5*wn_base, 0.5);
% zeta = 0.9;
% 
% sigma = zeta*wn;
% wd    = wn*sqrt(max(1-zeta^2,0));
% 
% if wd > 1e-8
%     p_dom = [-sigma + 1i*wd, -sigma - 1i*wd];
% else
%     p_dom = [-wn, -1.2*wn];
% end
% 
% p_extra = zeros(1,max(n-2,0));
% for k = 1:length(p_extra)
%     p_extra(k) = -(3+k)*wn;
% end
% 
% p_des = [p_dom, p_extra];
% Pdes  = real(poly(p_des));
% Pdes  = Pdes / Pdes(1);
% Pdes_sym = poly2sym(Pdes,s);
% 
% starts = [...
%      0.1   0.01;
%      1     0.1;
%      10    1;
%      100   10;
%     -0.1   0.01;
%     -1     0.1;
%     -10    1;
%      0.1  -0.01;
%      1    -0.1;
%      10   -1;
%     -0.1  -0.01;
%     -1    -0.1;
%     -10   -1];
% 
% bestJ   = inf;
% bestKp  = NaN;
% bestKd  = NaN;
% bestDen = [];
% 
% opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',500,'MaxFunEvals',3000);
% 
% for ii = 1:size(starts,1)
%     x0 = starts(ii,:);
% 
%     try
%         xopt = fminsearch(@(x)costo_pd(x,numGp,denGp,Pdes),x0,opts);
%         Kp_i = xopt(1);
%         Kd_i = xopt(2);
% 
%         den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp_i,Kd_i);
%         rts = roots(den_cl);
% 
%         if all(real(rts) < -1e-5)
%             J = costo_pd([Kp_i Kd_i],numGp,denGp,Pdes);
%             if J < bestJ
%                 bestJ   = J;
%                 bestKp  = Kp_i;
%                 bestKd  = Kd_i;
%                 bestDen = den_cl;
%             end
%         end
%     catch
%     end
% end
% 
% if isfinite(bestJ)
%     Kp_num = bestKp;
%     Kd_num = bestKd;
%     Plc_sym = poly2sym(bestDen,s);
% else
%     warning('No se encontró una solución estable para Kp y Kd.');
%     Kp_num = NaN;
%     Kd_num = NaN;
%     Plc_sym = sym(0);
% end
% 
% end
% 
% function J = costo_pd(x,numGp,denGp,Pdes)
% Kp = x(1);
% Kd = x(2);
% 
% if any(~isfinite([Kp Kd]))
%     J = 1e20;
%     return
% end
% 
% den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd);
% 
% if isempty(den_cl) || any(~isfinite(den_cl)) || abs(den_cl(1)) < 1e-12
%     J = 1e20;
%     return
% end
% 
% den_cl = real(den_cl / den_cl(1));
% 
% n = max(length(den_cl),length(Pdes));
% den_cl_pad = [zeros(1,n-length(den_cl)), den_cl];
% Pdes_pad   = [zeros(1,n-length(Pdes)),   Pdes];
% 
% err_coeff = norm(den_cl_pad - Pdes_pad,2)^2;
% 
% rts = roots(den_cl);
% pen_est = 1e6 * sum(max(real(rts)+1e-5,0).^2);
% pen_kd = 1e-3 / (abs(Kd) + 1e-4);
% pen_gain = 1e-8 * (Kp^2 + Kd^2);
% 
% J = err_coeff + pen_est + pen_kd + pen_gain;
% end
% 
% function den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd)
% Gc_num = [Kd Kp];
% 
% num_ol = conv(Gc_num,numGp);
% den_ol = denGp;
% 
% n = max(length(num_ol),length(den_ol));
% num_ol = [zeros(1,n-length(num_ol)), num_ol];
% den_ol = [zeros(1,n-length(den_ol)), den_ol];
% 
% den_cl = den_ol + num_ol;
% den_cl = quitar_ceros_inicio(real(den_cl));
% end
% 
% function p = quitar_ceros_inicio(p)
% idx = find(abs(p) > 1e-12,1,'first');
% if isempty(idx)
%     p = 0;
% else
%     p = p(idx:end);
% end
% end