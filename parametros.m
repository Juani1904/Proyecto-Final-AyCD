%SCRIPT de almacenamiento de parametros - Propiedades fisicas del sistema

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
kwu = 236e6;            % Rigidez unitaria a tracción del cable de izaje [N/m]
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
%% CN2 - Datos generales
T_s2=0.001; %Tiempo muestreo control nivel 2

%% CN2 - Controlador de movimiento - Carro
p1_t=0;
p2_t=-bt/Mt;
w_pos_t=4*p2_t;
n_t=2;
J_eqt = Jtm_ttb + (it^2 * Jtd) / rtd^2;
b_eqt = btm + (it^2 * btd) / rtd^2 + bt;
b_ta = (n_t * w_pos_t * rtd * J_eqt) / it - (rtd / it) * b_eqt;
K_tsa = (n_t * w_pos_t^2 * J_eqt * rtd) / it;
K_tsia = (w_pos_t^3 * J_eqt * rtd) / it;


%% CN2 - Controlador de movimiento - Izaje
J_eqh=(1/rhd)*(2*Jhd_hEb+2*Jhm_hb*(ih)^2);
b_eqh=(1/rhd)*(2*bhd+2*bhm*(ih)^2);
p1_h=0;
p2_h=-b_eqh/J_eqh;
w_pos_h=6*p2_h;
n_h=2;
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
