%SCRIPT de almacenamiento de parametros - Propiedades fisicas del sistema
syms b_ta K_tsa K_tsia b_ha K_hsa K_hsia
%-----------------------------------------------DATOS CONSIGNA----------------------------------------------------------------------
%% Datos generales
Yt0 = 45;               % Altura fija de las poleas de izaje en el carro [m]
Hc = 2.59;              % Altura del container estándar [m]
Wc = 2.44;              % Ancho del container estándar [m]
Ms = 15000;             % Masa del spreader + headblock (sin container) [kg]
Mc_max = 50000;         % Masa máxima del container totalmente cargado [kg]
Mc_min = 2000;          % Masa mínima del container vacío [kg]
g = 9.80665;            % Aceleración gravitatoria [m/s^2]
%% Carga apoyada - Parametros de contacto
Kcy = 1.8e9;            % Rigidez de compresión por contacto vertical [N/m]
bcy = 10.0e6;           % Fricción interna (amortiguamiento vertical) [N/(m/s)]
bcx = 1.0e6;            % Fricción de arrastre horizontal por contacto [N/(m/s)]
%% Cable de acero (wirerope) de IZAJE equivalente - Parametros unitarios
kwu = 236e6; h           % Rigidez unitaria a tracción del cable de izaje [N/m]
bwu = 150;              % Fricción interna unitaria a tracción del cable [N/(m/s)]
Lh0=110;                % Longitud de despliegue fijo del wirerope de izaje (Tambor-Extremo fijo), en adicion a las 2 partes colgantes variable lh(t)
%% Accionamiento del Sistema de Izaje
rhd = 0.75;                   % Radio primitivo del tambor (una sola corrida de cable) [m]
Jhd_hEb = 3800;               % Momento de inercia equivalente del eje lento [kg*m^2]
bhd = 8.0;                    % Coef. de fricción mecánica viscosa del eje lento [N*m/(rad/s)]
bhEb = 2.2e9;                 % Coef. de fricción viscosa del freno de emergencia cerrado [N*m/(rad/s)]
ThEb_Max = 1.1e6;             % Torque máximo del freno de emergencia cerrado [N*m]
ih = 22.0;                    % Relación de transmisión total de la caja reductora (i:h = 22:1)
Jhm_hb = 30.0;                % Momento de inercia equivalente del eje rápido [kg*m^2]
bhm = 18.0;                   % Coef. de fricción mecánica viscosa del eje rápido [N*m/(rad/s)]
bhb = 100e6;                  % Coef. de fricción viscosa del freno de operación cerrado [N*m/(rad/s)]
Thb_Max = 50.0e3;             % Torque máximo del freno de operación cerrado [N*m]
tauhm = 1.0e-3;               % Constante de tiempo del modulador de torque en el motor [s]
Thm_Max = 20.0e3;             % Torque máximo de motorización/frenado regenerativo del motor [N*m]
vh_max_cc = 1.5;
vh_max_sc = 3;
ah_max = 0.75;
M_x = 25000;

%% Carro y cable de acero (wirerope) de CARRO equivalente
Mt = 30000;                  % Masa equivalente del carro y componentes asociados [kg]
bt = 90.0;                   % Coef. de fricción mecánica viscosa equivalente del carro [N/(m/s)]
Ktw = 480e3;                 % Rigidez equivalente del cable tensado del carro [N/m]
btw = 3.0e3;                 % Fricción interna total del cable tensado del carro [N/(m/s)]

%% Accionamiento de Traslación del Carro
rtd = 0.50;                  % Radio primitivo del tambor (1 sola corrida de cable) [m]
Jtd = 1200;                  % Momento de inercia del eje lento (tambor y salida de reductora) [kg*m^2]
btd = 1.8;                   % Fricción mecánica viscosa del eje lento [N*m/(rad/s)]
it = 30.0;                   % Relación de transmisión total de la caja reductora (i:t = 30:1)
Jtm_ttb = 7.0;               % Momento de inercia del eje rápido (motor y entrada de reductora) [kg*m^2]
btm = 6.0;                   % Fricción mecánica viscosa del eje rápido [N*m/(rad/s)]
btb = 5.0e6;                 % Fricción viscosa del freno de operación cerrado [N*m/(rad/s)]
Ttb_Max = 5.0e3;             % Torque máximo del freno de operación cerrado [N*m]
tautm = 1.0e-3;              % Constante de tiempo del modulador de torque (motor del carro) [s]
Ttm_Max = 4.0e3;             % Torque máximo de motorización/frenado regenerativo del motor [N*m]
vt_max = 4;
at_max = 0.8;

%-----------------------------------------------DATOS AGREGADOS----------------------------------------------------------------------
%% Asignacion aleatoria de la masa del carro (Entorno)
%Mc_Xvect=
Mc_X=25000;
%% Perfil de obstaculos (Entorno)
N=ceil(80/Wc);
dmax=ceil(Yt0/2); %Altura maxima de apilado (Asumida)
step = Hc; %Paso de apilado
valores = 0:step:dmax; %Posibles valores de 0 a dmax, con paso step
y0 = valores(randi(numel(valores),1,N)); %Vector de perfil de obstaculos
y0(15)=0; %Setea en cero la altura correspondiente al borde entre el agua y el muelle
%% CN2 - Datos generales
T_s2=0.001; %Tiempo muestreo control nivel 2

%% CN2 - Controlador de movimiento - Carro
% p1_t=0;
% p2_t=bt/Mt;
% w_pos_t=500*p2_t;
% n_t=2.5; 
% % J_eqt = Jtm_ttb + (it^2 * Jtd) / rtd^2;
% % b_eqt = btm + (it^2 * btd) / rtd^2 + bt;
% J_eqt = Mt + Jtm_ttb * ((it^2) / rtd^2) + (Jtd/rtd^2)+Ms+M_x;
% b_eqt = btm *((it^2) / rtd^2)+  btd/(rtd^2) + bt;
% b_ta = ((n_t * w_pos_t * rtd * J_eqt) / it) - ((rtd / it) * b_eqt);
% %K_tsa = (n_t * w_pos_t^2 * J_eqt * rtd) / it;
% K_tsa = w_pos_t * b_ta;
% %K_tsia = (w_pos_t^3 * J_eqt * rtd) / it;
% K_tsia = w_pos_t * K_tsa;

J_eqt = Jtm_ttb*it^2/rtd^2 + Jtd/rtd^2 + Ms + M_x;   
b_eqt = btm*it^2/rtd^2 + btd/rtd^2;

coef_t = [rtd*(J_eqt+Mt)/it; rtd*(b_eqt+bt)/it; 0];
polos_t = roots(coef_t)

w_pos_t = 10*polos_t(2)
n_t = 2.9;

a = (rtd*bt/it + rtd*b_eqt/it + b_ta)/(rtd*(Mt+J_eqt)/it);
b = K_tsa/(rtd*(Mt+J_eqt)/it);
c = K_tsia/(rtd*(Mt+J_eqt)/it);
Gt_tranf = [a b c];
polin_sint_serie_t = [n_t*w_pos_t  n_t*w_pos_t^2  w_pos_t^3];

soluct = solve([Gt_tranf==polin_sint_serie_t], [b_ta, K_tsa, K_tsia]);
b_ta = -double(soluct.b_ta)
K_tsa = double(soluct.K_tsa)
K_tsia = -double(soluct.K_tsia)

roots([1; (rtd*bt/it + rtd*b_eqt/it + b_ta)/(rtd*(Mt+J_eqt)/it); K_tsa/(rtd*(Mt+J_eqt)/it); K_tsia/(rtd*(Mt+J_eqt)/it)])


%% CN2 - Controlador de movimiento - Izaje
% J_eqh=(1/rhd)*(2*Jhd_hEb+2*Jhm_hb*(ih)^2);
% b_eqh=(1/rhd)*(2*bhd+2*bhm*(ih)^2);
% p1_h=0;
% p2_h=-b_eqh/J_eqh;
% w_pos_h=6*p2_h;
% n_h=2;

J_eqh = (1/rhd)*(2*Jhd_hEb+2*Jhm_hb*(ih)^2);  
b_eqh = (1/rhd)*(2*bhd+2*bhm*(ih)^2);

coef_h = [-(1/ih)*(J_eqh+M_x*rhd/2); -b_eqh/ih; 0];
polos_h = roots(coef_h)

w_pos_h = 10*polos_h(2)
n_h = 2;

d = ih*(b_eqh/ih - b_ha)/(rhd*M_x/2+J_eqh);
e = -K_hsa*ih/(rhd*M_x/2+J_eqh);
f = -K_hsia*ih/(rhd*M_x+J_eqh);
Gh_tranf = [d e f];
polin_sint_serie_h = [n_h*w_pos_h  n_h*w_pos_h^2  w_pos_h^3];

soluch = solve([Gh_tranf==polin_sint_serie_h], [b_ha, K_hsa, K_hsia]);
b_ha = -double(soluch.b_ha)
K_hsa = double(soluch.K_hsa)
K_hsia = -double(soluch.K_hsia)

roots([1; ih*((b_eqh/ih) - b_ha)/((rhd*M_x/2)+J_eqh); -K_hsa*ih/((rhd*M_x/2)+J_eqh); -K_hsia*ih/(rhd*M_x+J_eqh)])

%% CN2 - Controlador de movimiento - Oscilacion Carga


%% CN2 -Modulador de Torque equivalente - Motor-Drive Izaje
A_hm=-1/tauhm;
B_hm= 1/tauhm;
C_hm= 1;
D_hm= 0;

%% CN2 -Modulador de Torque equivalente - Motor-Drive Carro
A_tm=-1/tautm;
B_tm= 1/tautm;
C_tm= 1;
D_tm= 0;


%% CN1 - Datos generales
T_s1=0.02; %Tiempo muestreo control nivel 2
dt=T_s1;
Fhw_min=1;
whm_min=10;
wtm_min=10;
xt_ref=-14; %ya n   o la toma de aca
dxt_max=4;
v_t_max=dxt_max;
v_t_max_2=dxt_max;
a_t_max=0.75;
xx_c0 =-30:0.01:50;
a_l_max=0.75;

%% Condiciones Iniciales
%En subs_acc_carro
dx_td_ini = 0; 
x_td_ini = 0;
%En subs_tras_carro
dx_t_ini=0;
x_t_ini=0;

%En subs_carga
dx_l_ini=0;
x_l_ini=0;
dy_l_ini=0;
y_l_ini=30;

%En acc_izaje
dl_h_ini=0;
l_h_ini=Yt0-y_l_ini;

% %% Linealizacion Jacobiana 
% 
% syms xt vt th w Ftw real
% syms Mt ml l bt g real
% 
% % Estados y entrada
% x = [xt; vt; th; w];
% u = Ftw;
% 
% % Definiciones auxiliares
% Dx  = Mt + ml*sin(th)^2;
% Nx = Ftw - bt*vt + ml*l*w^2*sin(th) + ml*g*sin(th)*cos(th);
% Nth = -(Mt+ml*sin(th)^2)*g*sin(th)-cos(th)*(Ftw-bt*vt+ml*l*w^2*sin(th)+ml*g*sin(th)*cos(th));
% Dth = l*(Mt+ml*sin(th)^2);
% 
% % Dinámica
% f1 = vt;
% f2 = Nx / Dx;
% f3 = w;
% f4 = Nth / Dth;
% 
% f = [f1; f2; f3; f4];
% 
% % Jacobianos
% A = jacobian(f, x);
% B = jacobian(f, u);
% 
% % (opcional) simplificación
% A = simplify(A);
% B = simplify(B);
% 
% % Modelo LPV evaluado en el punto de operacion
% 
% xt0 = 0;
% w0 = 0;
% Mt0 = M_x;
% ml0 = Mc_X;
% bt0 = btw;
% g0 = g;
% 
% A0 = subs(A, {xt,w}, {xt0,w0});
% B0 = subs(B, {xt,w}, {xt0,w0});
% 
% % Selección de salida y vector que representa la acción de delta-theta
% C = [0 1 0 0];            % salida y = delta v_t
% 
% %Para encontrar la funcion de transferencia
% SI_A0 = s*eye(4) - A0;
% Gp = simplify( C * inv(SI_A0)* B(:,1));  % TF de la planta
% 
% %% Controlador de carro
% syms Ktp Kti Ktd real
% 
% % Obtencion de Funcion de transferencia
% Gct = Ktp+Kti/s+Ktd*s;
% 
% %% Controlador de balanceo
% 
% % Obtencion de Funcion de transferencia
% syms Kp Kd s real
% 
% Gc = Kp+Kd*s;
% 
% Gcb = (Gct+Gc)*Gp/(1+(Gc+Gct)*Gp);
% Gcb = simplify(Gcb);
% [~, denGcb] = numden(Gcb);
% denGcb = expand(denGcb);
% %% Diseño por amortiguamiento crítico
% % Usamos una wn mayor a la de lazo abierto
% wn = 1.5*8;
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
% % % Forma robusta: extraer polinomio en vector
% % [~,P_des_d] = numden(P_des);
% % poly_des = expand(P_des_d);
% 
% % Igualación por coeficientes (grado 4 -> 5 ecuaciones)
% coeffs_lc = coeffs(P_lc,s,'All');
% coeffs_des = coeffs(P_des, s, 'All');
% 
% % % alinear longitudes rellenando con ceros si es necesario
% % nL = length(coeffs_lc); nD = length(coeffs_des);
% % n = max(nL, nD);
% % coeffs_lc = [zeros(1,n-nL), coeffs_lc];
% % coeffs_des = [zeros(1,n-nD), coeffs_des];
% 
% % formar sistema de ecuaciones (vectorial)
% eqs = coeffs_lc - coeffs_des;
% eqs = simplify(eqs);           % ecuaciones simbólicas
% 
% % Resolver simbólicamente Kp, Kd
% sol = solve(eqs, [Ktp Kti Ktd Kp Kd]);
% 
% % recoger soluciones posibles
% Kp_sol = sol.Kp;
% Kd_sol = sol.Kd;