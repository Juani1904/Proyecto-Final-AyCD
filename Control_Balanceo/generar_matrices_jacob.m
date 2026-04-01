function [Kp_mat, Kd_mat] = generar_matrices_jacob( ...
    xt0,vt0,th0,w0,Ftw0,Mt0,ml_vec,l_vec,bt0,g0)

nml = length(ml_vec);
nl  = length(l_vec);

Kp_mat = NaN(nml,nl);
Kd_mat = NaN(nml,nl);

for i = 1:nml
    for j = 1:nl
        ml0 = ml_vec(i);
        l0  = l_vec(j);

        try
            [Kp, Kd] = jacob_gs(xt0,vt0,th0,w0,Ftw0,Mt0,ml0,l0,bt0,g0);

            Kp_mat(i,j) = Kp;
            Kd_mat(i,j) = Kd;

            fprintf('ml=%.3f, l=%.3f -> Kp=%.6f, Kd=%.6f\n', ...
                ml0, l0, Kp, Kd);

        catch ME
            warning('Falló cálculo para ml = %.3f, l = %.3f. Motivo: %s', ...
                ml0, l0, ME.message);
        end
    end
end

end