% Reads a mesh from the specified file in the same fashion as Simple2D (FEM-Silvester).
function [numNodes, Nodes, numElements, Elements, numBC, BC] = ReadMesh(fileName)

% Open the file for reading only.
fid = fopen(fileName, 'r');

% 1. Read lines corresponding to node and coordinates.
% 1.1. Read the comment line and display it.
while 1
    tline = fgetl(fid);
    if ~isequal(tline,''),   
        disp('-->> Reading:')
        disp(tline)
        break;
    else
        continue;
    end
end

% 1.2. Read the number of NODES to process.
numNodes = fscanf(fid,'%d',1);
% 1.3. Read the NODES given by numNodes.
%Nodes = textscan(fid, '%n %n', numNodes);
Nodes = fscanf(fid, '%f %f', [2 numNodes]);
Nodes = Nodes';

% 2. Read lines corresponding to ELEMENTS and distributed sources.
% 2.1. Read the comment line and display it.
while 1
    tline = fgetl(fid);
    if ~isequal(tline,''),   
        disp('-->> Reading:')
        disp(tline)
        break;
    else
        continue;
    end
end

% tline = fgetl(fid);
disp('-->> Reading:')
disp(tline)
% 2.2. Read the number of ELEMENTS to process.
numElements = fscanf(fid,'%d',1);
% 2.3. Read the numbers of ELEMENTS given by N.
% Elements = textscan(fid, '%d %d %d %n', numElements);
Elements = fscanf(fid, '%d %d %d %f', [4 numElements]);
Elements = Elements';

% 3. Read lines corresponding to BOUNDARY CONDITIONS.
% 3.1. Read the comment line and display it.
while 1
    tline = fgetl(fid);
    if ~isequal(tline,''),   
        disp('-->> Reading:')
        disp(tline)
        break;
    else
        continue;
    end
end

% 3.2. Read the number of BOUNDARY CONDITIONS to process.
numBC = fscanf(fid,'%d',1);
% 3.3. Read the numbers of BOUNDARY CONDITIONS given by N.
% BC = textscan(fid, '%d %n', numBC);
BC = fscanf(fid, '%d %f', [2 numBC]);
BC = BC';

% Close file
fclose(fid);
