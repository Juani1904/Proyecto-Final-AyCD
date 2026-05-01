%Añado los paths de las funciones
addpath("Control_Balanceo\")
%SCRIPT de almacenamiento de parametros -  Propiedades fisicas del sistema
syms b_ta K_tsa K_tsia b_ha K_hsa K_hsia
%-----------------------------------------------DATOS CONSIGNA----------------------------------------------------------------------
%% Datos generales
Yt0 = 45;               % Altura fija de las poleas de izaje en el carro [m]
Hc = 2.59;              % Altura del container estándar [m]
Wc = 2.44;              % Ancho del container estándar [m]
Ms = 15000;             % Masa del spreader + headblock (sin container) [kg]
Mc_max = 50000;         % Masa máxima del container totalmente cargado [kg]
Mc_min = 2000;          % -----Masa mínima del container vacío [kg]
g = 9.80665;            % Aceleración gravitatoria [m/s^2]
a_t_max=0.8;            % Aceleración máxima del carro
a_h_max=0.75;           % Aceleración maxima de izaje

%% Carga apoyada - Parametros de contacto
Kcy = 1.8e9;            % Rigidez de compresión por contacto vertical [N/m]
bcy = 10.0e6;           % Fricción interna (amortiguamiento vertical) [N/(m/s)]
bcx = 1.0e6;            % Fricción de arrastre horizontal por contacto [N/(m/s)]
%% Cable de acero (wirerope) de IZAJE equivalente - Parametros unitarios
kwu = 236e6;           % Rigidez unitaria a tracción del cable de izaje [N/m]
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
Jtm_tb = 7.0;               % Momento de inercia del eje rápido (motor y entrada de reductora) [kg*m^2]
btm = 6.0;                   % Fricción mecánica viscosa del eje rápido [N*m/(rad/s)]
btb = 5.0e6;                 % Fricción viscosa del freno de operación cerrado [N*m/(rad/s)]
Ttb_Max = 5.0e3;             % Torque máximo del freno de operación cerrado [N*m]
tautm = 1.0e-3;              % Constante de tiempo del modulador de torque (motor del carro) [s]
Ttm_Max = 4.0e3;             % Torque máximo de motorización/frenado regenerativo del motor [N*m]
vt_max = 4;
at_max = 0.8;

%-----------------------------------------------DATOS AGREGADOS----------------------------------------------------------------------
%% Perfil de obstaculos (Entorno)
N=ceil(80/Wc);
dmax=ceil(Yt0/4); %Altura maxima de apilado (Asumida)
step = Hc; %Paso de apilado
valores = 0:step:dmax; %Posibles valores de 0 a dmax, con paso step
y0 = valores(randi(numel(valores),1,N)); %Vector de perfil de obstaculos
y0(15)=5; %Setea en cero la altura correspondiente al borde entre el agua y el muelle

% Generacion de perfil de obstaculos para apilamiento de containers
% Muelle a la izquierda (contenores aislados, llegada por camiones)
% Barco a la derecha (apilamiento por bloques / mesetas)

% N = ceil(80/Wc);
% dmax = ceil(Yt0/2);
% step = Hc;
% 
% % Indice del borde agua–muelle
% idx_borde = ceil(30/Wc);
% 
% % Definicion de zonas
% quay_idx = 1:(idx_borde-1);      % muelle (izquierda)
% ship_idx = (idx_borde+1):N;      % barco (derecha)
% 
% % Alturas maximas por zona (en niveles discretos)
% maxShip = max(2, floor(0.90*dmax/step));
% maxQuay = max(1, floor(0.60*dmax/step));
% 
% % Parametros de correlacion espacial
% Lship = 4;   % bloques mas largos en barco
% Lquay = 2;   % menor correlacion en muelle
% 
% % Señal correlacionada simple
% corrSignal = @(n,L) filter(ones(1,max(1,L))/max(1,L), 1, randn(1,n));
% 
% % Normalizacion suave a [0,1] para evitar extremos
% to01 = @(z) max(0, min(1, 0.5 + 0.20*z));
% 
% % ----- BARCO: apilamiento por bloques -----
% zShip = corrSignal(numel(ship_idx), Lship);
% uShip = to01(zShip);
% 
% gammaShip = 1.6;    % controla sesgo hacia alturas medias
% lvlShip = round((uShip.^(1/gammaShip)) * maxShip);
% 
% % Suavizado adicional para generar mesetas
% lvlShip = round(filter([1 2 3 2 1]/9, 1, lvlShip));
% lvlShip = max(0, min(maxShip, lvlShip));
% 
% % ----- MUELLE: contenedores aislados -----
% zQuay = corrSignal(numel(quay_idx), Lquay);
% uQuay = to01(zQuay);
% 
% gammaQuay = 2.4;    % mayor sesgo a alturas bajas
% lvlQuay = round((uQuay.^(1/gammaQuay)) * maxQuay);
% lvlQuay = max(0, min(maxQuay, lvlQuay));
% 
% % Pequeña irregularidad vertical
% jitterProb = 0.10;
% jitter = (rand(size(lvlQuay)) < jitterProb) .* (randi(3, size(lvlQuay)) - 2);
% lvlQuay = max(0, min(maxQuay, lvlQuay + jitter));
% 
% % Mascara de ocupacion para generar stacks aislados
% p_single = 0.12;    % stacks individuales
% p_pair   = 0.03;    % pares ocasionales
% 
% nq = numel(lvlQuay);
% mask = false(1, nq);
% 
% i = 1;
% while i <= nq
%     r = rand;
%     if r < p_pair && i < nq
%         mask(i) = true;
%         mask(i+1) = true;
%         i = i + 3;      % dejo un hueco para evitar bloques largos
%     elseif r < (p_pair + p_single)
%         mask(i) = true;
%         i = i + 2;      % dejo un hueco
%     else
%         i = i + 1;
%     end
% end
% 
% lvlQuay(~mask) = 0;
% 
% % Evito triples consecutivos dejando solo el central
% for k = 2:(nq-1)
%     if lvlQuay(k)>0 && lvlQuay(k-1)>0 && lvlQuay(k+1)>0
%         lvlQuay(k-1) = 0;
%         lvlQuay(k+1) = 0;
%     end
% end
% 
% % ----- Construccion del perfil final -----
% y0 = zeros(1, N);
% y0(quay_idx) = lvlQuay * step;
% y0(ship_idx) = lvlShip * step;
% 
% % Altura cero en el borde agua–muelle
% y0(idx_borde) = 0;
% 
% % Limite al salto maximo entre columnas para evitar picos no realistas
% maxStepJump = 2*step;
% for k = 2:N
%     if y0(k) > y0(k-1) + maxStepJump
%         y0(k) = y0(k-1) + maxStepJump;
%     elseif y0(k) < y0(k-1) - maxStepJump
%         y0(k) = y0(k-1) - maxStepJump;
%     end
% end
% 
% y0(idx_borde) = 0;

%% Asignacion aleatoria de la masa del carro (Entorno)
Mc_X=randi([Mc_min,Mc_max]); %Este valor mas adelante va a desaparecer
Mc_Xvect=randi([Mc_min, Mc_max], 1, N); %Generación de masas contenedores aleatoria entre Mc min y Mc max
Mc_Xvect_sobrecarga=ones(1,N)*(Mc_max+1000);

%% CN2 - CONTROL REGULATORIO

% ------------------------------------------------CN2 - DATOS GENERALES---------------------------------------------------

T_s2=0.001; %Tiempo muestreo control nivel 2

% ---------------------------------------CN2 - CONTROLADOR DE MOVIMIENTO CARRO -----------------------------------------
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

J_eqt = Jtm_tb*it^2/rtd^2 + Jtd/rtd^2 + Ms + M_x;   
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


% ------------------------------------------------CN2 - CONTROLADOR DE MOVIMIENTO IZAJE -------------------------------
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


% --------------------------------CN2 - MODULADOR DE TORQUE EQUIVALENTE - MOTOR DRIVE IZAJE -------------------------
A_hm=-1/tauhm;
B_hm= 1/tauhm;
C_hm= 1;
D_hm= 0;

% --------------------------------CN2 - MODULADOR DE TORQUE EQUIVALENTE - MOTOR DRIVE IZAJE -----------------------
A_tm=-1/tautm;
B_tm= 1/tautm;
C_tm= 1;
D_tm= 0;

%% --------------------------------------------CN2 - CONTROL BALANCEO -----------------------------------------------

%Condiciones iniciales
xt0  = -20;
%vt0  = 0;
%th0  = 0;
w0   = 0;
%Ftw0 = 0;


ml_vec = linspace(Ms+Mc_min,Ms+Mc_max,8);
ml_vec = [Ms,ml_vec]; %Se explicita que el primer valor de Ml debe ser el escenario para spreader vacio, masa Ms
l_vec  = linspace(5,45,10);
vt_vec = linspace(0,vt_max,10);
at_vec = linspace(-at_max,at_max,10);

% GS = preparar_gain_scheduling_lookup(xt0,vt_vec,at_vec,w0,Mt,ml_vec,l_vec,bt,g,true);


%[Kp,Kd,A0,B0,Gp_sym,Pdes_sym,Plc_sym] = jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,2000,10,bt0,g0);


%% CN1 - CONTROL SUPERVISOR

%--------------------------------------------------CN1 - DATOS GENERALES---------------------------------------------

%Tiempo muestreo control nivel 1
T_s1=0.02;
F01_Ts1=T_s1;
F02_Ts1=T_s1;
F03_Ts1=T_s1;

%Fuerza de izaje min (Determinacion cable flojo/tenso)
Fhw_min=10;

%Velocidad angular de izaje minima
whm_min=10;

%Velocidad angular de carro minima
wtm_min=10;

%COMPLETAR
xt_ref=-20; %ya no la toma de aca

%Aceleración máxima carro
F01_a_t_max=a_t_max;
F03_a_t_max=a_t_max;
%Completar
xx_c0 =-30:0.01:50;

%Velocidad máxima carro
dxt_max=4;
F03_dxt_max=dxt_max;

%Aceleración máxima izaje
F02_a_h_max=a_h_max;
F03_a_h_max=a_h_max;

%Angulo de izaje max tolerable, carro detenido
tita_l_max_estable=1*(pi/180); %1 grado (En radianes)

%% CN0 - CONTROL SEGURIDAD

%--------------------------------------------------CN0 - DATOS GENERALES---------------------------------------------

T_s0=0.02;

%% Condiciones Iniciales
%En subs_acc_carro
dx_td_ini = 0;      %MODIFICABLE
x_td_ini = -20;     %MODIFICABLE

%En subs_tras_carro
dx_t_ini=dx_td_ini;
x_t_ini=x_td_ini;

%En subs_carga
dx_l_ini=0;         %MODIFICABLE
x_l_ini=x_td_ini;   
dy_l_ini=0;         %MODIFICABLE
y_l_ini=20;         %MODIFICABLE

%En acc_izaje
dl_h_ini=dy_l_ini;
l_h_ini=Yt0-y_l_ini;







% %% Linealizacion Jacobiana 
% 
% syms s xt vt th w Ftw real
% syms Mt_sym ml_sym l_sym bt_sym g_sym real
% 
% % Estados y entrada
% x = [xt; vt; th; w];
% u = Ftw;
% 
% % Definiciones auxiliares
% Dx  = Mt_sym + ml_sym*sin(th)^2;
% Nx = Ftw - bt_sym*vt + ml_sym*l_sym*w^2*sin(th) + ml_sym*g_sym*sin(th)*cos(th);
% Nth = -(Mt_sym+ml_sym*sin(th)^2)*g_sym*sin(th)-cos(th)*(Ftw-bt_sym*vt+ml_sym*l_sym*w^2*sin(th)+ml_sym*g_sym*sin(th)*cos(th));
% Dth = l_sym*(Mt_sym+ml_sym*sin(th)^2);
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
% A0 = subs(A, {xt,w,Mt_sym,ml_sym,bt_sym,g_sym}, {xt0,w0,Mt0,ml0,bt0,g0});
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
% syms Kp Kd real
% 
% Gc = Kp+Kd*s;
% 
% Gcb = (Gc+Gct)*Gp/(1+(Gc+Gct)*Gp);
% Gcb = simplify(Gcb);
% [~, denGcb] = numden(Gcb);
% denGcb = simplify(denGcb);
% 
% %% Diseño por amortiguamiento crítico
% % % Usamos una wn mayor a la de lazo abierto
% wn = 1.5*8;
% w_pos_cb = 10*wn;
% 
% % Polinomio dominante de segundo orden crítico
% P_dom = s^2 + 2*wn*s + wn^2;
% 
% % Polos no dominantes (rápidos)
% alpha = 5*wn;
% beta = 6*wn;
% P_des = expand(P_dom * (s + alpha) * (s + beta));
% 
% % Coeficientes del polinomio del lazo cerrado
% coeffs_gcb=coeffs(denGcb, s, 'All');
% 
% % Coeficientes del polinomio del lazo cerrado
% coeffs_des=coeffs(P_dom, s, 'All');
% 
% % % Igualación de coeficientes
% % Pdiff = (expand(D_lc - D_des));
% % coeffs_vec = coeffs(Pdiff, s);
% % eqs = coeffs_vec == 0;
% 
% % % % Forma robusta: extraer polinomio en vector
% % % [~,P_des_d] = numden(P_des);
% % % poly_des = expand(P_des_d);
% 
% % % Igualación por coeficientes (grado 4 -> 5 ecuaciones)
% % coeffs_lc = coeffs(P_lc,s,'All');
% % coeffs_des = coeffs(P_des, s, 'All');
% 
% % alinear longitudes rellenando con ceros si es necesario
% nL = length(coeffs_gcb); nD = length(coeffs_des);
% n = max(nL, nD);
% coeffs_gcb = [zeros(1,n-nL), coeffs_gcb];
% coeffs_des = [zeros(1,n-nD), coeffs_des];
% % 
% % % formar sistema de ecuaciones (vectorial)
% % eqs = coeffs_gcb - coeffs_des;
% % eqs = simplify(eqs);           % ecuaciones simbólicas
% 
% Gcb_tranf = coeffs_gcb;
% 
% polin_sint_serie_cb = coeffs_des;
% 
% soluch = solve([Gcb_tranf==polin_sint_serie_cb], [Kp, Kd]);
% 
% Kp = double(soluch.Kp);
% Kd = double(soluch.Kd);
% 
% % % Resolver simbólicamente Kp, Kd
% % sol = solve(eqs==0, [Kp Kd]);
% % 
% % % recoger soluciones posibles
% % Kp_sol = sol.Kp;
% % Kd_sol = sol.Kd;