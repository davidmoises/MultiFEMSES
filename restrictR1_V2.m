function R2h = restrictR1_V2(R, Nodes, Nodes2h, h_parts)
% Restricts the R from a fine mesh (Nodes) to a coarser mesh (Nodes2h)
% Fine partition
limF = h_parts*3+1;

% Coarse partition
limC = round((h_parts/2)*3+1);

% Restric R vales to coarse mesh
R2h = Nodes2h;

%% Restrict
for i=1:limC
    for j=1:limC
        indexC = (i-1)*limC + j;
        if R2h(indexC,1) == 0,
           indexF = (i-1)*2*limF + (j*2 - 1);
           R2h(indexC,2) = (1/16)*(4*R(indexF) + 2*(R(indexF - limF)...
               + R(indexF + limF) + R(indexF + 1) + R(indexF - 1)) + ...
               (R(indexF - limF - 1) + R(indexF - limF + 1) + ...
                R(indexF + limF - 1) + R(indexF + limF + 1)));
        end
    end
end