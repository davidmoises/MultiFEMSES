function [e2h, rp2h, iter] = smootherC_V2(r2h, h, x_part, tol, max_iter)
% Smoother for a given coarse mesh with uniform step size h.

iter = 0;
err = inf;
s = size(r2h);
e2h = zeros(s);
e2h(:,1) = r2h(:,1);
rp2h = zeros(s);
rp2h(:,1) = r2h(:,1);
% Fine partition
lim = x_part*3+1;

%% Smoother cicle for the error.
while (err>tol && iter<=max_iter)
    iter = iter + 1;
    % Store previos error value.
    e2h_ant = e2h(:,2);
    % Detemernime new error value.
    for i=1:s(1),
        if e2h(i,1)==0,
           e2h(i,2) = 0.25*(e2h(i + 1,2) + e2h(i - 1,2) + e2h(i + lim,2)...
                             + e2h(i - lim,2) - r2h(i,2)*(h^2));
        end
    end
    % Compute error on the current iterate.
    if(iter==1)
       err_n=norm(e2h(:,2)); 
    end
    err = norm(e2h(:,2)-e2h_ant)/err_n;
end

%% Compute the residual of the current solution per node for the nest
% coarser grid.
for i=1:s(1),
    if rp2h(i,1)==0,
%         rp2h(i,2) = r2h(i,2) - (1/h^2)*(e2h(i + 1,2) + e2h(i - 1,2)...
%                     + e2h(i + lim,2) + e2h(i - lim,2) - 4*e2h(i,2));
        rp2h(i,2) = -(1/h^2)*(e2h(i + 1,2) + e2h(i - 1,2)...
                    + e2h(i + lim,2) + e2h(i - lim,2) - 4*e2h(i,2));
    end
end