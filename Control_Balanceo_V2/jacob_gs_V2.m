function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym, Gcb_sym] = ...
    jacob_gs_V2(xt0,vt0,at0,w0,Mt0,ml0,l0,bt0,g0)
warning('off', 'symbolic:solve:SolutionsDependOnConditions');

th0 = atan(-at0/g0);
Ftw0 = at0*ml0*sin(th0)^2 - l0*ml0*sin(th0)*w0^2 + Mt0*at0 - (g0*ml0*sin(2*th0))/2;
z=1; %Parametro Zeta (Amortiguacion Critica)

%% Parámetros matriz A

a22 = -bt0/(ml0*sin(th0)^2 + Mt0);

a23 = (ml0*(l0*cos(th0)*w0^2 + g0*cos(2*th0)))/(ml0*sin(th0)^2 + Mt0) - (2*ml0*cos(th0)*sin(th0)*(Ftw0 + ml0*(l0*sin(th0)*w0^2 + (g0*sin(2*th0))/2) - bt0*vt0))/(ml0*sin(th0)^2 + Mt0)^2;

a24 = (2*l0*ml0*sin(th0)*w0)/(ml0*sin(th0)^2 + Mt0);

a42 = (bt0*cos(th0))/(l0*(ml0*sin(th0)^2 + Mt0));

a43 = (sin(th0)*(Ftw0 + ml0*(l0*sin(th0)*w0^2 + (g0*sin(2*th0))/2) - bt0*vt0))/(l0*(ml0*sin(th0)^2 + Mt0)) - (g0*cos(th0))/l0 - (ml0*cos(th0)*(l0*cos(th0)*w0^2 + g0*cos(2*th0)))/(l0*(ml0*sin(th0)^2 + Mt0)) + (2*ml0*cos(th0)^2*sin(th0)*(Ftw0 + ml0*(l0*sin(th0)*w0^2 + (g0*sin(2*th0))/2) - bt0*vt0))/(l0*(ml0*sin(th0)^2 + Mt0)^2);

a44 = -(2*ml0*cos(th0)*sin(th0)*w0)/(ml0*sin(th0)^2 + Mt0);

%% Parámetros matriz B

b21 = 1/(ml0*sin(th0)^2 + Mt0);

b41 = -cos(th0)/(l0*(ml0*sin(th0)^2 + Mt0));

%% Matrices de estado

A = [0    1   0   0;
     0   a22 a23 a24;
     0    0   0   1;
     0   a42 a43 a44];

B = [ 0;
     b21;
      0;
     b41];
%% Calculo de Kp y Kd
rad = (A(2,3)*B(4) - A(4,3)*B(2)) / B(2);
wn = sqrt(rad);     %frecuencia natural
Kd_num = -(A(4,3)*B(2)*B(4) - A(2,3)*B(4)^2 - A(2,2)*A(2,4)*B(4)^2 - A(4,2)*A(4,4)*B(2)^2 + B(2)*B(4)*wn^2 - 2*A(4,2)*B(2)^2*wn*z + A(2,2)*A(4,4)*B(2)*B(4) + A(2,4)*A(4,2)*B(2)*B(4) + 2*A(2,2)*B(2)*B(4)*wn*z) ...
        /(A(2,2)^2*B(4)^2 - 2*A(2,2)*A(4,2)*B(2)*B(4) + 2*z*A(2,2)*B(4)^2*wn + A(4,2)^2*B(2)^2 - 2*z*A(4,2)*B(2)*B(4)*wn + B(4)^2*wn^2);
den = (B(2) + Kd_num*B(4));
Kp_num = -(((A(2,4)*B(4) - A(4,4)*B(2) - Kd_num*A(2,2)*B(4) + Kd_num*A(4,2)*B(2))/den - 2*wn*z) * (B(2) + Kd_num*B(4))) / B(4);



%Ftw0 = (Mt0 + ml0)*at0 + bt0*vt0;

% syms xt vt th w Ftw Mt ml l bt g dd_xt dd_th s real
% 
% %% Estados e input
% x = [xt; vt; th; w]; %Estados x1=xt x2=d_xt=vt x3=th x4=d_th=w
% u = Ftw;
% 
% %% Modelo no lineal exacto - Modelo dinámico
% eq1 = (Mt + ml)*dd_xt + ml*l*dd_th*cos(th) - ml*l*w^2*sin(th) == Ftw - bt*vt;
% eq2 = l*dd_th + dd_xt*cos(th) + g*sin(th) == 0;
% 
% sol = solve([eq1, eq2], [dd_xt, dd_th]);
% 
% dd_xt_expr  = simplify(sol.dd_xt);
% dd_th_expr = simplify(sol.dd_th);
% 
% %% Dinamica de estados: xdot = f(x,u) MODELO EN ESPACIO DE ESTADOS
% f = [vt; dd_xt_expr; w; dd_th_expr];
% 
% %% Jacobianos
% A = simplify(jacobian(f, x));
% B = simplify(jacobian(f, u));
% 
% %% Punto de operacion
% solF = solve(dd_xt_expr,Ftw);
% Ftw0= subs(solF,[vt; th; w ; Mt; ml; l; bt; g],[vt0; th0; w0; Mt0; ml0; l0; bt0; g0]);
% 
% vars_sym = [xt; vt; th; w; Ftw; Mt; ml; l; bt; g];
% vals_num = [xt0; vt0; th0; w0; Ftw0; Mt0; ml0; l0; bt0; g0];
% 
% A0 = subs(A, vars_sym, vals_num);
% B0 = subs(B, vars_sym, vals_num);
% 
% A0_num = double(A0);
% B0_num = double(B0);
% 
% %% Planta local linealizada Gp(s) = Theta(s)/Vt(s)
% Cth = [0 0 1 0];
% Cv  = [0 1 0 0];
% 
% GthF = minreal(tf(ss(A0_num, B0_num, Cth, 0)));
% GvF  = minreal(tf(ss(A0_num, B0_num, Cv, 0)));
% 
% Gp_tf = minreal(GthF / GvF);
% 
% [numGp, denGp] = tfdata(Gp_tf, 'v');
% numGp = quitar_ceros_inicio(real(numGp));
% denGp = quitar_ceros_inicio(real(denGp));
% 
% if numel(denGp) ~= 3 || numel(numGp) ~= 2
%     error('La sintonia serie implementada requiere Gp(s) de orden 2 con numerador de orden 1.');
% end
% 
% Gp_sym = poly2sym(numGp, s) / poly2sym(denGp, s);
% 
% %% Coeficientes de Gp(s) = (b1*s + b0)/(a2*s^2 + a1*s + a0)
% a2 = denGp(1);
% a1 = denGp(2);
% a0 = denGp(3);
% 
% b1 = numGp(1);
% b0 = numGp(2);
% 
% %% Parametros de sintonia serie sobre Gcb(s)
% % Gcb(s) = Gp(s) / (1 + Gc(s)*Gp(s))
% % Gc(s)  = Kd*s + Kp
% %
% % Denominador de Gcb:
% % Dcb(s) = Dp(s) + Gc(s)Np(s)
% %        = (a2 + Kd*b1)s^2 + (a1 + Kd*b0 + Kp*b1)s + (a0 + Kp*b0)
% %
% % Se impone el polinomio deseado normalizado:
% % Pdes(s) = s^2 + alpha*s + beta
% % con alpha = n_p*w_pos_p y beta = w_pos_p^2
% 
% polos_gp = roots(denGp);
% polos_gp = polos_gp(isfinite(polos_gp));
% 
% if isempty(polos_gp)
%     w_nat_gp = sqrt(abs(a0/a2));
% else
%     w_nat_gp = max(abs(polos_gp));
% end
% 
% lambda_ss = 1.3;
% w_pos_p   = lambda_ss * w_nat_gp;
% n_p       = 2.8;
% 
% alpha = n_p * w_pos_p;
% beta  = w_pos_p^2;
% 
% %% Calculo de Kp y Kd igualando el denominador de Gcb al polinomio deseado
% % Ecuaciones:
% % (a1 + Kd*b0 + Kp*b1)/(a2 + Kd*b1) = alpha
% % (a0 + Kp*b0)/(a2 + Kd*b1)         = beta
% 
% M = [b1,            (b0 - alpha*b1);
%      b0,            (-beta*b1)];
% 
% rhs = [alpha*a2 - a1;
%        beta*a2  - a0];
% 
% if abs(det(M)) < 1e-12 || rcond(M) < 1e-12
%     error('El sistema para calcular Kp y Kd esta mal condicionado.');
% end
% 
% solK = M \ rhs;
% 
% Kp_num = real(solK(1));
% Kd_num = real(solK(2));
% 
% if ~isfinite(Kp_num) || ~isfinite(Kd_num)
%     error('No se pudieron obtener ganancias finitas.');
% end
% 
% %% Controlador PD
% Gc_num = [Kd_num Kp_num];
% 
% %% Funcion de transferencia a lazo cerrado
% num_cb = numGp;
% den_cb = suma_poly(denGp, conv(Gc_num, numGp));
% 
% if isempty(den_cb) || abs(den_cb(1)) < 1e-12
%     error('El denominador de lazo cerrado es invalido.');
% end
% 
% den_cb = real(den_cb / den_cb(1));
% num_cb = real(num_cb / den_cb(1)); %#ok<NASGU>
% 
% Plc_sym  = poly2sym(den_cb, s);
% Gcb_sym  = poly2sym(numGp, s) / poly2sym(den_cb, s);
% 
% %% Polinomio deseado
% Pdes = [1 alpha beta];
% Pdes_sym = poly2sym(Pdes, s);
% 
% %% Chequeo de estabilidad
% r_cl = roots(den_cb);
% if any(real(r_cl) >= 0)
%     warning('El lazo cerrado calculado no quedo estable para este punto de operacion.');
% end
% 
% end
% 
% function p = quitar_ceros_inicio(p)
% tol = 1e-12;
% idx = find(abs(p) > tol, 1, 'first');
% if isempty(idx)
%     p = 0;
% else
%     p = p(idx:end);
% end
% end
% 
% function p = suma_poly(p1, p2)
% n = max(length(p1), length(p2));
% p1 = [zeros(1, n-length(p1)), p1];
% p2 = [zeros(1, n-length(p2)), p2];
% p  = quitar_ceros_inicio(real(p1 + p2));
% end


