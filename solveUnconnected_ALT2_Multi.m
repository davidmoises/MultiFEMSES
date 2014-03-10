% Solucionador de equaciones utilizando bloques independientes de elementos
% y estilo Jacobi de solucion. Se utilizan los elementos fuera de la
% diagonal para actualizar cada iteracion.
function [x, iter, err_norm] = solveUnconnected_ALT2_Multi(numElements, Sdis_J, b_dis, Cst, tolerance, max_iter, x_guess)
x = x_guess;
iter = 0;
err_norm = inf;
temp_x = x_guess;

while (err_norm >= tolerance)
    x_old = x;
    temp_x = zeros(numElements*3,1);
    bcount = 1;

    for i=1:numElements
%         if Sdis_J{i,5}~=0,
%             % Get the corresponding x-vector elements.
%             x_local = x_old(Sdis_J{i,3});    
%             % Solve a new iteration and average with the previous one.
%             temp = b_dis{i} - (Sdis_J{i,2})*x_local; % Sdis_J{1,2} corresponds to the off diagonal elements.
%             temp_x = [temp_x ; Sdis_J{i,1}*temp];  % Sdis_J{1,1} corresponds to the diagonal elements only ALREADY INVERTED.
%         end
try
        % Get the corresponding x-vector elements.
        x_local = x_old(Sdis_J{i,3});
        % Solve a new iteration and average with the previous one.
%         temp = b_dis(i) - Sdis_J{i,2}*x_local; % Sdis_J{1,2} corresponds to the off diagonal elements.
        temp_x(bcount:bcount+2) = Sdis_J{i,1}*(b_dis{i} - Sdis_J{i,2}*x_local);
catch
%    [i, x_local, b_dis(bcount:bcount+2,2), Sdis_J{i,2}];
   i, x_local, b_dis{i}, Sdis_J{i,2}
end
        bcount = bcount + 3;
    end

    % Connect the 2 solutions.
    x = Cst*temp_x;
%     iter_i = size(Cst,1);
%     for i=1:iter_i
%         x(i) = 0;
%         iter_j=length(Cst{i,1});
%         for j=1:iter_j
%             x(i) = x(i) + temp_x(Cst{i,1}(j))*Cst{i,2}(j);
%         end
%     end
 
    % Check for convergence
    iter = iter + 1;
try
    err_norm = norm(x_old-x);
%     err_norm = norm(x_old-temp_x); % NEW on 02262014
    if (iter>max_iter)
       break; 
    end
catch
    iter, length(x_old), length(temp_x)
    error('ERROR: in function solveUnconnected_ALT2_Multi');
end

end
