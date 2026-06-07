% Posicion de referencia
xt_ref_tiempo = out.xt_ref.Time;
xt_ref_amplitud = out.xt_ref.Data;
xt_ref_data = [xt_ref_tiempo, xt_ref_amplitud];

% Posicion real
x_t_real_tiempo = out.x_t_real.Time;
x_t_real_amplitud = out.x_t_real.Data;
x_t_real_data = [x_t_real_tiempo, x_t_real_amplitud];

% Longitud de referencia
lh_ref_tiempo = out.lh_ref.Time;
lh_ref_amplitud = out.lh_ref.Data;
lh_ref_data = [lh_ref_tiempo, lh_ref_amplitud];

% Longitud real
l_h_real_tiempo = out.l_h_real.Time;
l_h_real_amplitud = out.l_h_real.Data;
l_h_real_data = [l_h_real_tiempo, l_h_real_amplitud];

% Velocidad de referencia de carro
dxt_ref_tiempo = out.dxt_ref.Time;
dxt_ref_amplitud = out.dxt_ref.Data;
dxt_ref_data = [dxt_ref_tiempo, dxt_ref_amplitud];

% Velocidad real de carro
dx_t_real_tiempo = out.dx_t_real.Time;
dx_t_real_amplitud = out.dx_t_real.Data;
dx_t_real_data = [dx_t_real_tiempo, dx_t_real_amplitud];

% Velocidad de referencia de cable
dlh_ref_tiempo = out.dlh_ref.Time;
dlh_ref_amplitud = out.dlh_ref.Data;
dlh_ref_data = [dlh_ref_tiempo, dlh_ref_amplitud];

% Velocidad real de cable
dl_h_real_tiempo = out.dlh_ref.Time;
dl_h_real_amplitud = out.dlh_ref.Data;
dl_h_real_data = [dl_h_real_tiempo, dl_h_real_amplitud];

% Aceleracion de referencia
ddxt_ref_tiempo = out.ddxt_ref.Time;
ddxt_ref_amplitud = out.ddxt_ref.Data;
ddxt_ref_data = [ddxt_ref_tiempo, ddxt_ref_amplitud];

% Angulo de referencia
tita_ref_tiempo = out.tita_ref.Time;
tita_ref_amplitud = out.tita_ref.Data;
tita_ref_data = [tita_ref_tiempo, tita_ref_amplitud];

% Angulo real
tita_real_tiempo = out.tita_real.Time;
tita_real_amplitud = out.tita_real.Data;
tita_real_data = [tita_real_tiempo, tita_real_amplitud];


writematrix(xt_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\xt_ref_simulacion.csv');
writematrix(x_t_real_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\x_t_real_simulacion.csv');
writematrix(lh_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\lh_ref_simulacion.csv');
writematrix(l_h_real_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\l_h_real_simulacion.csv');
writematrix(dxt_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\dxt_ref_simulacion.csv');
writematrix(dx_t_real_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\dx_t_real_simulacion.csv');
writematrix(dlh_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\dlh_ref_simulacion.csv');
writematrix(dl_h_real_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\dl_h_real_simulacion.csv');
writematrix(ddxt_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\ddxt_ref_simulacion.csv');
writematrix(tita_ref_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\tita_ref_simulacion.csv');
writematrix(tita_real_data, 'C:\Users\juani\OneDrive\Escritorio\Analisis graficas proyectos\CSV\tita_real_simulacion.csv');