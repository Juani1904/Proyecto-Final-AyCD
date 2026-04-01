function Mout = rellenar_nan_matriz(Min)

Mout = Min;

if all(isnan(Mout),'all')
    error('La matriz completa es NaN. No se puede interpolar.');
end

[nr,nc] = size(Mout);

[X,Y] = meshgrid(1:nc,1:nr);

idx_valid = ~isnan(Mout);
idx_nan   = isnan(Mout);

% Interpolación lineal sobre los puntos válidos
Mout(idx_nan) = griddata( ...
    X(idx_valid), Y(idx_valid), Min(idx_valid), ...
    X(idx_nan),   Y(idx_nan), ...
    'linear');

% Si quedaron NaN en bordes/regiones externas, completar con nearest
idx_nan2 = isnan(Mout);
if any(idx_nan2,'all')
    Mout(idx_nan2) = griddata( ...
        X(idx_valid), Y(idx_valid), Min(idx_valid), ...
        X(idx_nan2),  Y(idx_nan2), ...
        'nearest');
end

end

