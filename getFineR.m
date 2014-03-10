function r = getFineR(numElements, Sdis, b_dis, Cs, x)
% Compute the residual for the fine mesh.

count = 0;
temp_r = zeros(numElements*3,1);
for i=1:numElements
    % Get the corresponding x-vector elements.
    x_local = x(Sdis{i,2});
%     A_local = Sdis{i,1};

    % Solve a new iteration and average with the previous one.
%     r_local = b_dis{i} - A_local*x_local ;
    temp_r(count+1:count+3) = b_dis{i} - Sdis{i,1}*x_local ;
%     temp_r = [temp_r; r_local];
    count = count + 3;
end

% Connect the 2 solutions.
r = Cs*temp_r;
%     iter_i = size(Cst,1);
%     for i=1:iter_i
%         x(i) = 0;
%         iter_j=length(Cst{i,1});
%         for j=1:iter_j
%             x(i) = x(i) + temp_r(Cst{i,1}(j))*Cst{i,2}(j);
%         end
%     end
