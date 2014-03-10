% Functions to create the Disconnected Coeficient matrix for 2D 1sr order 
% triangular systems and the Conection matrix.
% Returns a reduced system according to the BC defined.
function [Sdis, C, CNodes, b_dis] = CreateDiscon3D_NoEmptyEl_V3_Multigrid(numElements, Elements, numNodes, Nodes, numBC, BC)
% Define some local variables.
count = 0;
x=zeros(3,1);
y=zeros(3,1);
% Stores: a) the reduced Se matrix, b) the nodes that are not BC in the 
% order of apperance in the element, and c) the size(# of non BC nodes).
% e.i. Sdis{1,:} = {Reduced Se matrix, elementNodes-list, num-non-BC-nodes}
% <- for element 1.
Sdis = cell(numElements,3);
% C: Conection matrix : [Se (element-nodes) (# of nodes in Se)]
C = sparse(numElements*3, numNodes);
% CNodes: stores the row and column indices in the Connection matrix using
% the mesh node numbering.
CNodes = cell(1,2);% Added on 4-2-2010
CNodes{1,1} = zeros(numElements*3,1);% Numbering of the Rows in C. Added on 4-2-2010
CNodes{1,2} = zeros(numNodes,1);% Numbering of the Columns in C. Added on 4-2-2010
b_dis = cell(numElements,1);
disNodeCount = 0;

% Iterate over ELEMENTS to build local S matrix
for i=1:numElements
    localBCnodes = zeros(3,1);
    b_BC_temp = [];
    tempElements = (Elements(i,(1:3)))';% Added on 4-2-2010
    % Iterate over the Element nodes to determine their coordinates and
    % connetivity.
    nodeBCcount = 0;
    bval = zeros(3,1);
    for j=1:3
        % Set the coordinate values for each of the three nodes in the
        % triangula element.
        tempNode = Elements(i,j);
        x(j) = Nodes(tempNode,1);
        y(j) = Nodes(tempNode,2);
        
        % Determine if current node belongs to a BC.
        for k=1:numBC
            tempBC = BC(k,1);
            if (tempNode == tempBC)
                localBCnodes(j) = j;
                % Build local BC-RHS vector.
                bval(j) = BC(k,2);
                b_BC_temp = [b_BC_temp ; bval(j)];
                nodeBCcount = nodeBCcount + 1;
                break;
            end
        end

        % Set element connectivity in the C matrix
        C(count+j,Elements(i,j)) = 1;
        CNodes{1,1}(count+j) = Elements(i,j);% Added on 4-2-2010
        CNodes{1,2}(Elements(i,j)) = Elements(i,j);% Added on 4-2-2010
    end
    
    % Increase the count on elements that only contain BC nodes and thus
    % need to be removed.
    Se = zeros(3,3);
    if(nodeBCcount==3)
%         tempElements = -1*ones(3,1);
        Se = eye(3);
    else
        disNodeCount = disNodeCount + 3 - nodeBCcount;
        
        % Compute 4xTriangle area.
        % Note: if the triangle area is not positive then the result may change,
        % thus it is made to be positive always.
        t_area = 2*abs((x(2)*y(3) + x(1)*y(2) + x(3)*y(1)) - (x(1)*y(3) + x(2)*y(1) + x(3)*y(2)));

        % Build Se from Sij factors.        
        Se(1,1) = (y(2) - y(3))^2 + (x(3) - x(2))^2;
        Se(2,2) = (y(3) - y(1))^2 + (x(1) - x(3))^2;
        Se(3,3) = (y(1) - y(2))^2 + (x(2) - x(1))^2;
        Se(1,2) = (y(2) - y(3))*(y(3) - y(1)) + (x(3) - x(2))*(x(1) - x(3));
        Se(2,1) = Se(1,2);
        Se(1,3) = (y(2) - y(3))*(y(1) - y(2)) + (x(3) - x(2))*(x(2) - x(1));
        Se(3,1) = Se(1,3);
        Se(2,3) = (y(3) - y(1))*(y(1) - y(2)) + (x(1) - x(3))*(x(2) - x(1));
        Se(3,2) = Se(2,3);
        Se = Se./t_area;
        
        if (nodeBCcount~=0)
            % Build local vector to remove the independent rows of the Se matrix to
            % build.
            localBCnodes(localBCnodes == 0) = [];

            % Remove the independent rows from the matrix Se using localBCnodes.
            Se(localBCnodes,:) = zeros(length(localBCnodes),3); %Se(localBCnodes,:) = [];
            % Assign independent columns of Se to the RHS matrix b_dis
            bval = bval-Se(:,localBCnodes)*b_BC_temp;
            % added on 3/4/2014
            Se(:,localBCnodes) = zeros(3,length(localBCnodes));
            Se(localBCnodes,localBCnodes) = eye(length(localBCnodes));
        end
    end
    b_dis{i} = bval;
    Sdis{i,1} = Se;
    Sdis{i,2} = tempElements;% Node numbers of the current element. Added on 4-2-2010
    Sdis{i,3} = localBCnodes; %length(tempElements);% # of nodes added into Sdis. Added on 4-2-2010

    % Increase count.
    count = count + 3;
end