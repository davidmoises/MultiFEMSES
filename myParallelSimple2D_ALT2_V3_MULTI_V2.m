% myParallelSimple2D_ALT2_V3_MULTI_V2
% Author: David Fernandez
% Universidad del Zulia
% My second Multigrid implementation
clear all; clc;

%% Variable declaration
% fileName = './MeshGen/MsquareCoax2_small.msh'; x_part = 2;
fileName = './MeshGen/MsquareCoax4_small.msh'; x_part = 4; fileName2h= './MeshGen/MsquareCoax2_small.msh';
% fileName = './MeshGen/MsquareCoax8_small.msh'; x_part = 8; fileName2h= './MeshGen/MsquareCoax4_small.msh';
% fileName = './MeshGen/MsquareCoax16_small.msh'; x_part = 16; fileName2h= './MeshGen/MsquareCoax8_small.msh';
% fileName = './MeshGen/MsquareCoax32_small.msh'; x_part = 32; fileName2h= './MeshGen/MsquareCoax16_small.msh';
% fileName = './MeshGen/MsquareCoax64_small.msh'; x_part = 64; fileName2h= './MeshGen/MsquareCoax32_small.msh';
% fileName = './MeshGen/MsquareCoax128_small.msh'; x_part = 128; fileName2h= './MeshGen/MsquareCoax64_small.msh';
y_part = x_part;
x_bounds = [0 0.02/3 0.02*2/3 0.02];
y_bounds = x_bounds;
h = (x_bounds(2)-x_bounds(1))/x_part;

%% Read Nodes, Elements and Boundary Conditions from file.
tic;
[numNodes, Nodes, numElements, Elements, numBC, BC] = ReadMesh(fileName);
[numNodes2, Nodes2, numElements2, Elements2, numBC2, BC2] = ReadMesh(fileName2h);% Remove
[numNodesh, Nodesh] = ReadNodesM(fileName);
[numNodes2h, Nodes2h] = ReadNodesM(fileName2h);

%% Clear data if not needed.
% clear Sdis Nodes Elements BC C CNodes
% disp('Finished Clearing data!!!')

%% Solve unconnected system.
tic;
tolerance = 1e-5;
max_iter = 10000;
x_guess = zeros(numNodes,1);
no_multigrid = 0;
use_vcycle = 1;

if no_multigrid,
    %% Compute the Disconnected Coeficient Sdis matrix and the Connection matrix.
    [Sdis, C, CNodes, b_dis] = CreateDiscon3D_NoEmptyEl_V2_Multigrid(numElements, ...
        Elements, numNodes, Nodes, numBC, BC); % This is the one working for BuildFEM2DSysJ_Multi y MultiV2
    
    %% Build new J-system to solve.
    [Cs, Sdis_J] = BuildFEM2DSysJ_Multi(C, Sdis, numElements); % This one works
    Cs = Cs';
    Assembly_time = toc;

    %% Basic FEMSES unconnected solver.
    [x, iter, err_norm] = solveUnconnected_ALT2_Multi(numElements, Sdis_J, b_dis, ...
                                        Cs, tolerance, max_iter, x_guess);
    disp(['Iterations for FEMSES(solveUnconnected_ALT2): ' num2str(iter)]);
    fprintf('With err-norm of: %f \n', err_norm);
    
    %% Display vector results.
    % Show results.
    % disp('x-->>')
    % disp(x)
    idx = find(x==0);
    x(idx)=[];
%     [sol2D] = MyPlotMesh(10,0, x, x_bounds, y_bounds, x_part, y_part);
%     [C h] = contour(sol2D);
else
    %% Compute the Disconnected Coeficient Sdis matrix and the Connection matrix.
%     [Sdis, C, CNodes, b_dis] = CreateDiscon3D_NoEmptyEl_V2_Multigrid(numElements, ...
%         Elements, numNodes, Nodes, numBC, BC);
    [Sdis, C, CNodes, b_dis] = CreateDiscon3D_NoEmptyEl_V3_Multigrid(numElements, ...
        Elements, numNodes, Nodes, numBC, BC);

    %% Build new J-system to solve.
%     [Cs, Sdis_J] = BuildFEM2DSysJ_Multi(C, Sdis, numElements);
    [Cs, Sdis_J] = BuildFEM2DSysJ_MultiV2(C, Sdis, numElements);
    Cs = Cs';
    Assembly_time = toc;
    
    if use_vcycle==1,
        % Basic V-cycle
        max_cycle_iter = 50;
        % Solve the fine graind with FEMSES relaxation for a tolerance of
        % 0.1.
        [xf, iter1, err_norm] = solveUnconnected_ALT2_Multi(numElements, ... 
                                Sdis_J, b_dis, Cs, 0.01, 100, x_guess); %tolerance*1000
        fprintf('Vcycle: iter->%i err-nom->%f \n', iter1, err_norm);
        % Restrict the  solution to the coarser grid (2-levels only).
        r = getFineR(numElements, Sdis, b_dis, Cs, xf);
        r2h = restrictR1_V2(r,Nodes,Nodes2h, x_part);
        % Smoother at the coarse grid
        [e2h, rp2h, iter2] = smootherC_V2(r2h, h*2, x_part/2, tolerance, max_cycle_iter);
        % Interpolate back to the fine grid
        ef = interpolateI2(e2h, Nodesh, x_part);
        x_guess = correctMapV2(ef, xf);
        % Final smoother on the fine grid.
        [x, iter3, err_norm] = solveUnconnected_ALT2_Multi(numElements, ...
                                Sdis_J, b_dis, Cs, tolerance, max_iter, x_guess); 
        disp('Iterations for MULTI-FEMSES(solveUnconnected_ALT2): ');
        iterF = iter1 + iter3
        iterC = iter2
        [iter1 iter2 iter3]
        
        idx = find(abs(x-0)<tolerance);
        x(idx)=[];
        idx = find(abs(x-10)<tolerance);
        x(idx)=[];
    else
        % Use Full Multigrid
        max_cycle_iter = 50;
        r2h = Nodes2h;
        [x2h, rp2h, iter0] = smootherC_FMG(r2h, Nodes2h, h*2, x_part/2, tolerance, max_cycle_iter);
        xc = interpolateI1(x2h, Nodesh, numNodesh, x_part);
        x_guess = correctMap(xc, Sdis_J, Nodesh, numNodesh, x_guess);
        [xf, iter1, err_norm] = solveUnconnected_ALT2(numElements, Sdis_J, b_dis, ...
                                            Cs, tolerance*1000, 100, x_guess);
        r = getFineR(newNumElements, Sdis_J, b_dis, Cs, xf);
        r2h = restrictR1(r,Nodes,Nodes2h, x_part);
        [e2h, rp2h, iter2] = smootherC(r2h, Nodes2h, h*2, x_part/2, tolerance, max_cycle_iter);
        ef = interpolateI1(e2h, Nodesh, numNodesh, x_part);
        x_guess = correctMap(ef, Sdis_J, Nodesh, numNodesh, xf);
        [x, iter3, err_norm] = solveUnconnected_ALT2(numElements, Sdis_J, b_dis, ...
                                            Cs, tolerance, max_iter, x_guess); 
        disp('Iterations for FULL-MULTI-FEMSES(solveUnconnected_ALT2): ');
        iterF = iter1 + iter3
        iterC = iter2 + iter0
        [iter0 iter1 iter2 iter3]
    end
end
New_Solver_time = toc;

%% Display Triangular Mesh
% Elements = Elements(:,1:3);
% z_val = zeros(size(Nodes, 1), 1);
% trisurf(Elements,Nodes(:,1),Nodes(:,2),z_val);

%% Display Equipotential plot of Results
figure;
[sol2D] = MyPlotMesh(10,0, x, x_bounds, y_bounds, y_part, x_part);
[C, h] = contour(sol2D);
disp('done!!!')
return
%% Determine errors
x1 = x;
% load('matrix2');
load('matrix4');
% load('matrix8');
% load('matrix16');
% load('matrix32');
% load('matrix64');
% load('matrix128');

err = norm(x1 - x);
w_new = 0.5*(x1'*(A*x1));
w_old = 0.5*(x'*(A*x));
w_err = abs((w_new-w_old)/w_old)*100;