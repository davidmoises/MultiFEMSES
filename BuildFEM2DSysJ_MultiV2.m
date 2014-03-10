%% Build new J-system to solve.
function [Cs, Sdis_J] = BuildFEM2DSysJ_MultiV2(C, Sdis, numElements)
Sdis_J = cell(numElements,4);
% Identify the scaling factors for each element matrix to scale the
% Connection matrix.
iter_i = size(C,2);

% Compute the accunulated sum of all diagonal matrix entries corresponding
% to superposed nodes identified as repeated ones in the columns of C.
CScaling = zeros(iter_i,2); % [Accum; count]
for i=1:numElements
    for j=1:3,
        val = Sdis{i,1}(j,j); 
        index = Sdis{i,2}(j); 
        CScaling(index,1) = CScaling(index,1) + val;
        CScaling(index,2) = CScaling(index,2) + 1;
    end
end

% Now copy these scaling factors into the C matrix
Cs = C;
[index_i, index_j] = find(C);
N = length(index_i);
for i=1:N % each element in C
    % Scale the new "Cs" matrix elements depending on the scaling factor
    % matrix "CScaling"
    if (CScaling(index_j(i),2)>=2 && CScaling(index_j(i),1)~=0)
        Cs(index_i(i),index_j(i)) = 1/CScaling(index_j(i),1);
    else
        Cs(index_i(i),index_j(i)) = 1/CScaling(index_j(i),2);
    end
end

% Separate the main diagonal elements from the diagonal ones to build Sdis_J.
count = 0;
for i=1:numElements
    Sdis_J{i,1} = diag(diag(Sdis{i,1}));
    Sdis_J{i,2} = Sdis{i,1} - Sdis_J{i,1};
    % Get the inverse of the main diagonal elements only.
%     Sdis_J{i,1} = inv(Sdis_J{i,1}); % deleted because of 0's in the diagonal
    for j=1:3,
        val = Sdis{i,1}(j,j);
        if (val~=0),
            index = Sdis{i,2}(j); 
            Cs(count+j,index) = Cs(count+j,index) * val;
            Sdis_J{i,1}(j,j) = 1/Sdis_J{i,1}(j,j);
        end
        Sdis_J{i,4}(j) = count+j; % global distributed nodes % MOVED ON 3/3/2014
    end
    Sdis_J{i,3} = Sdis{i,2}; % global node indices for this element
    count = count + 3;
end
