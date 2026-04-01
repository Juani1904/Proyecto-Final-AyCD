function GS = preparar_gain_scheduling_lookup( ...
    xt0,vt0,th0,w0,Ftw0,Mt0,ml_vec,l_vec,bt0,g0,aplicar_suavizado)

if nargin < 11
    aplicar_suavizado = false;
end

%% 1) Generar matrices crudas
[Kp_raw, Kd_raw] = generar_matrices_jacob( ...
    xt0,vt0,th0,w0,Ftw0,Mt0,ml_vec,l_vec,bt0,g0);

%% 2) Rellenar NaN
Kp_fill = rellenar_nan_matriz(Kp_raw);
Kd_fill = rellenar_nan_matriz(Kd_raw);

%% 3) Suavizado opcional
if aplicar_suavizado
    Kp_use = suavizar_matriz(Kp_fill);
    Kd_use = suavizar_matriz(Kd_fill);
else
    Kp_use = Kp_fill;
    Kd_use = Kd_fill;
end

%% 4) Guardar en estructura
GS = struct();
GS.ml_vec  = ml_vec(:);
GS.l_vec   = l_vec(:);
GS.Kp_raw  = Kp_raw;
GS.Kd_raw  = Kd_raw;
GS.Kp_fill = Kp_fill;
GS.Kd_fill = Kd_fill;
GS.Kp      = Kp_use;
GS.Kd      = Kd_use;

%% 5) Variables listas para Simulink
Kp_table = GS.Kp;
Kd_table = GS.Kd;
ml_breakpoints = GS.ml_vec;
l_breakpoints  = GS.l_vec;

assignin('base','GS',GS);
assignin('base','Kp_table',Kp_table);
assignin('base','Kd_table',Kd_table);
assignin('base','ml_breakpoints',ml_breakpoints);
assignin('base','l_breakpoints',l_breakpoints);

end

