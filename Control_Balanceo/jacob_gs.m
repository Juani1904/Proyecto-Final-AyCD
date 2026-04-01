function [Kp_num, Kd_num, A0_num, B0_num, Gp_sym, Pdes_sym, Plc_sym] = ...
    jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0)
% Ftw0 se mantiene solo por compatibilidad

syms xt vt th w Tm_a Tm real
syms Mt ml l bt g real
syms s real
syms tau_a K_a R_eq eta_eq real

x = [xt; vt; th; w; Tm_a];
u = Tm;

Tm0   = 0;
Tm_a0 = 0;

f_a = (-Tm_a + K_a*Tm)/tau_a;
Ftw = eta_eq*Tm_a/R_eq;

Dx  = Mt + ml*sin(th)^2;
Nx  = Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th);

Nth = -(Mt + ml*sin(th)^2)*g*sin(th) ...
      - cos(th)*(Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th));

Dth = l*(Mt + ml*sin(th)^2);

f1 = vt;
f2 = Nx / Dx;
f3 = w;
f4 = Nth / Dth;

f = [f1; f2; f3; f4; f_a];

A = simplify(jacobian(f,x));
B = simplify(jacobian(f,u));

vars_sym = {xt,vt,th,w,Tm_a,Tm,Mt,ml,l,bt,g,tau_a,K_a,R_eq,eta_eq};
vals_num = {xt0,vt0,th0,w0,Tm_a0,Tm0,Mt0,ml0,l0,bt0,g0,0.05,1,0.5,1};

A0 = simplify(subs(A, vars_sym, vals_num));
B0 = simplify(subs(B, vars_sym, vals_num));

A0_num = double(A0);
B0_num = double(B0);

C0 = [0 0 1 0 0];
D0 = 0;

Gp_ss = ss(A0_num,B0_num,C0,D0);
Gp_tf = minreal(tf(Gp_ss));

[numGp, denGp] = tfdata(Gp_tf,'v');
numGp = quitar_ceros_inicio(real(numGp));
denGp = quitar_ceros_inicio(real(denGp));

Gp_sym = poly2sym(numGp,s) / poly2sym(denGp,s);

n = length(denGp) - 1;

if n < 2
    warning('La planta local tiene orden menor que 2.');
    Kp_num = NaN;
    Kd_num = NaN;
    Pdes_sym = sym(0);
    Plc_sym  = sym(0);
    return
end

lambda = eig(A0_num);
real_part = abs(real(lambda));
real_part = real_part(real_part > 1e-4);

if isempty(real_part)
    wn_base = 0.5;
else
    wn_base = min(real_part);
end

wn   = max(1.5*wn_base, 0.5);
zeta = 0.9;

sigma = zeta*wn;
wd    = wn*sqrt(max(1-zeta^2,0));

if wd > 1e-8
    p_dom = [-sigma + 1i*wd, -sigma - 1i*wd];
else
    p_dom = [-wn, -1.2*wn];
end

p_extra = zeros(1,max(n-2,0));
for k = 1:length(p_extra)
    p_extra(k) = -(3+k)*wn;
end

p_des = [p_dom, p_extra];
Pdes  = real(poly(p_des));
Pdes  = Pdes / Pdes(1);
Pdes_sym = poly2sym(Pdes,s);

starts = [...
     0.1   0.01;
     1     0.1;
     10    1;
     100   10;
    -0.1   0.01;
    -1     0.1;
    -10    1;
     0.1  -0.01;
     1    -0.1;
     10   -1;
    -0.1  -0.01;
    -1    -0.1;
    -10   -1];

bestJ   = inf;
bestKp  = NaN;
bestKd  = NaN;
bestDen = [];

opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',500,'MaxFunEvals',3000);

for ii = 1:size(starts,1)
    x0 = starts(ii,:);

    try
        xopt = fminsearch(@(x)costo_pd(x,numGp,denGp,Pdes),x0,opts);
        Kp_i = xopt(1);
        Kd_i = xopt(2);

        den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp_i,Kd_i);
        rts = roots(den_cl);

        if all(real(rts) < -1e-5)
            J = costo_pd([Kp_i Kd_i],numGp,denGp,Pdes);
            if J < bestJ
                bestJ   = J;
                bestKp  = Kp_i;
                bestKd  = Kd_i;
                bestDen = den_cl;
            end
        end
    catch
    end
end

if isfinite(bestJ)
    Kp_num = bestKp;
    Kd_num = bestKd;
    Plc_sym = poly2sym(bestDen,s);
else
    warning('No se encontró una solución estable para Kp y Kd.');
    Kp_num = NaN;
    Kd_num = NaN;
    Plc_sym = sym(0);
end

end

function J = costo_pd(x,numGp,denGp,Pdes)
Kp = x(1);
Kd = x(2);

if any(~isfinite([Kp Kd]))
    J = 1e20;
    return
end

den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd);

if isempty(den_cl) || any(~isfinite(den_cl)) || abs(den_cl(1)) < 1e-12
    J = 1e20;
    return
end

den_cl = real(den_cl / den_cl(1));

n = max(length(den_cl),length(Pdes));
den_cl_pad = [zeros(1,n-length(den_cl)), den_cl];
Pdes_pad   = [zeros(1,n-length(Pdes)),   Pdes];

err_coeff = norm(den_cl_pad - Pdes_pad,2)^2;

rts = roots(den_cl);
pen_est = 1e6 * sum(max(real(rts)+1e-5,0).^2);
pen_kd = 1e-3 / (abs(Kd) + 1e-4);
pen_gain = 1e-8 * (Kp^2 + Kd^2);

J = err_coeff + pen_est + pen_kd + pen_gain;
end

function den_cl = den_lazo_cerrado_pd(numGp,denGp,Kp,Kd)
Gc_num = [Kd Kp];

num_ol = conv(Gc_num,numGp);
den_ol = denGp;

n = max(length(num_ol),length(den_ol));
num_ol = [zeros(1,n-length(num_ol)), num_ol];
den_ol = [zeros(1,n-length(den_ol)), den_ol];

den_cl = den_ol + num_ol;
den_cl = quitar_ceros_inicio(real(den_cl));
end

function p = quitar_ceros_inicio(p)
idx = find(abs(p) > 1e-12,1,'first');
if isempty(idx)
    p = 0;
else
    p = p(idx:end);
end
end