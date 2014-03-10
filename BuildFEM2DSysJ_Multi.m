%% Build new J-system to solve.
function [Cs, Sdis_J] = BuildFEM2DSysJ_Multi(C, Sdis, numElements)
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
        if (val~=0),
            index = Sdis{i,2}(j); 
            CScaling(index,1) = CScaling(index,1) + val;
            CScaling(index,2) = CScaling(index,2) + 1;
        end
    end
end

% Now copy these scaling factors into the C matrix
Cs = C;
[iter_i, iter_j] = size(C);
for i=1:iter_i % Row iteration
    for j=1:iter_j % Column iteration
        % Scale the new "Cs" matrix elements depending on the scaling factor
        % matrix "CScaling"
        if (Cs(i,j)~=0)
            if (CScaling(j,2)>=2)
                Cs(i,j) = 1/CScaling(j,1);
            end
        end
    end
end

% Multiply the non 0's or 1's factors in the matrix Cs matrix by the 
% corresponding diagonal factor S(j,j) in the Sdis matrix.
count = 0;
% for i=1:numElements
%     for j=1:3,
%         % Scaling of the Cs non 1 elements according to the Sdis & Cs elements.
%         for k=1:3,  % new implementation multigrid
%             if (Cs(count+k,Sdis{i,2}(j))~=0 && Cs(count+k,Sdis{i,2}(j))~=1)
%                 Cs(count+k,Sdis{i,2}(j)) = Cs(count+k,Sdis{i,2}(j)) * Sdis{i,1}(j,j);
%                 Sdis_J{i,3}(j) = count+k; % global distributed nodes
%                 break;
%             end
%         end
%     end
%     count = count + 3;
% end

% commented ON 3/3/2014
% for i=1:numElements
%     for j=1:3,
%         val = Sdis{i,1}(j,j); 
%         if (val~=0),
%             index = Sdis{i,2}(j); 
%             Cs(count+j,index) = Cs(count+j,index) * val;
% %             Sdis_J{i,3}(j) = count+j; % global distributed nodes
%         end
%         Sdis_J{i,3}(j) = count+j; % global distributed nodes % MOVED ON 3/3/2014
%     end
%     count = count + 3;
% end

% Separate the main diagonal elements from the diagonal ones to build Sdis_J.
for i=1:numElements
    Sdis_J{i,1} = diag(diag(Sdis{i,1}));
    Sdis_J{i,2} = Sdis{i,1} - Sdis_J{i,1};
    % Get the inverse of the main diagonal elements only.
%     Sdis_J{i,1} = inv(Sdis_J{i,1}); 
    for j=1:3,
        val = Sdis{i,1}(j,j);
        if (val~=0),
            index = Sdis{i,2}(j); 
            Cs(count+j,index) = Cs(count+j,index) * val;
%             Sdis_J{i,3}(j) = count+j; % global distributed nodes
            Sdis_J{i,1}(j,j) = 1/Sdis_J{i,1}(j,j);
        end
        Sdis_J{i,4}(j) = count+j; % global distributed nodes % MOVED ON 3/3/2014
        % Removed ON 3/3/2014
%         if Sdis_J{i,1}(j,j)~=0,
%            Sdis_J{i,1}(j,j) = 1/Sdis_J{i,1}(j,j);
%         end
    end
    Sdis_J{i,3} = Sdis{i,2}; % global node indices for this element
    count = count + 3;
end
