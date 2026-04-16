function [Kp_mat, Kd_mat] = generar_matrices_jacob( ...
    xt0,vt_vec,at_vec,w0,Mt0,ml_vec,l_vec,bt0,g0)

nml = length(ml_vec);
nl  = length(l_vec);
nvt = length(vt_vec);
nat = length(at_vec);

Kp_mat = NaN(nml,nl);
Kd_mat = NaN(nml,nl);

%for i = 1:nml
    %ml0 = ml_vec(i);
    ml0=15000;
    for j = 1:nl
        l0  = l_vec(j);
        %for k = 1:nvt
            %vt0 = vt_vec(k);
            vt0 = 0;
            for l = 1:nat
                at0 = at_vec(l);
                try
                    [Kp, Kd] = jacob_gs(xt0,vt0,at0,w0,Mt0,ml0,l0,bt0,g0);
        
                    % Kp_mat(i,j,k,l) = Kp;
                    % Kd_mat(i,j,k,l) = Kd;
                     Kp_mat(j,l) = Kp;
                     Kd_mat(j,l) = Kd;
        
                    fprintf('ml=%.3f, l=%.3f, vt=%.3f, at=%.3f -> Kp=%.6f, Kd=%.6f\n', ...
                        ml0, l0,vt0,at0, Kp, Kd);
        
                catch ME
                    warning('Falló cálculo para ml = %.3f, l = %.3f, vt=%.3f, at=%.3f. Motivo: %s', ...
                        ml0, l0, vt0,at0, ME.message);
                end
                
            end
        %end
        
        
    end
%end

end