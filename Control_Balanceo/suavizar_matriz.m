function Ms = suavizar_matriz(M)

kernel = [1 2 1;
          2 4 2;
          1 2 1] / 16;

Ms = conv2(M, kernel, 'same');

% Mantener bordes originales para no distorsionar tanto
Ms(1,:)   = M(1,:);
Ms(end,:) = M(end,:);
Ms(:,1)   = M(:,1);
Ms(:,end) = M(:,end);

end