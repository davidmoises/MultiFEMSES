function ef = interpolateI2(e2h, Nodesh, h_part)
% Fine partition
limF = h_part*3+1;
% Coarse partition
limC = round((h_part/2)*3+1);
% Bese for Interpolated values 
ef = Nodesh;

%% Copy aligned values and perform horizontal/vertical/diagonal sweeps.
for i=1:(limC-1)
    for j=1:(limC-1)
        indexC = (i-1)*limC + j;
        indexF = (i-1)*2*limF + (j*2 - 1);
        % Copy overlapping values and do all sweeps.
        if e2h(indexC,1) == 0,
           ef(indexF,2) = e2h(indexC,2);
           ef(indexF + 1,2) = 0.5*(e2h(indexC,2) + e2h(indexC + 1,2));
           ef(indexF + limF,2) = 0.5*(e2h(indexC,2) + e2h(indexC + limC,2));
           ef(indexF + limF + 1,2) = 0.25*(e2h(indexC,2) + e2h(indexC + 1,2) ...
                 + e2h(indexC + limC,2) + e2h(indexC + limC + 1,2));
        else 
           % Horizontal sweep.
           if ef(indexF + 1,1) == 0 
              ef(indexF + 1,2) = 0.5*(e2h(indexC,2) + e2h(indexC + 1,2));
           end
           % Vertical sweep.
           if ef(indexF + limF,1) == 0 
              ef(indexF + limF,2) = 0.5*(e2h(indexC,2) + e2h(indexC + limC,2));
           end
           % Diagonal sweep.
           if ef(indexF + limF + 1,1) == 0 
              ef(indexF + limF + 1,2) = 0.25*(e2h(indexC,2) + e2h(indexC + 1,2) ...
                    + e2h(indexC + limC,2) + e2h(indexC + limC + 1,2));
           end 
        end
    end
end