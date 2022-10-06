%%% Cleaning memory %%% 
clear;
clc;
%%% 3D truss structure %%%

% --- Nodal coordinates --- %
%       xi yi zi 
nodes = [0 0 0;                         % node 1 
        0 1 0;                          % node 2    
        sqrt(3)/2 0.5 0;                % node 3 
        0.5/sqrt(3) 0.5 sqrt(2/3)];     % node 4 
n = size(nodes,1);                      % number of nodes

% --- Elements --- %
%           i  j  EA      we assume a nominal EAo=1 
% --- Elements --- %
%      i j EA adT   
elem = [1 4 1 0;          % element 1
        4 2 1 0.05;          % element 2
        4 3 1 0];         % element 3
m = size(elem,1);       % number of elements 

% --- Stiffness matrix --- %
K = zeros(3*n,3*n);     % initialization;  
fT = zeros(3*n,1);   % initialization 

% assembly; 
for q = 1:m 
    i = elem(q,1); xi = nodes(i,1); yi = nodes(i,2); zi = nodes(i,3); % extracting the coordinates of the nodes of pertinence 
    j = elem(q,2); xj = nodes(j,1); yj = nodes(j,2); zj = nodes(j,3); 
    L = norm([xj;yj;zj]-[xi;yi;zi]);  Lv(q)=L; 
    kq = elem(q,3)/L;  % stiffness 
    kv(q) = kq; 
    qv = ([xj;yj;zj]-[xi;yi;zi])/L; 
    Ne = (kq*[-qv;qv]); 
    Nev(q,:) = Ne; 
    Ke = 1/kq*Ne*Ne'; 
    edofs = [3*i-2,3*i-1,3*i,3*j-2,3*j-1,3*j];
    fT(edofs) = fT(edofs)+elem(q,4)*L*Ne;  % the vector fT relaed to element q is fT=aDTq*L*Ne, with aDTq =elem(q,4) (*) 
    K(edofs,edofs) = K(edofs,edofs)+Ke;  
end

% --- Boundary conditions --- %
fix_nod = [1 2 3];  % nodes that are fixed 
fix_dofs = [];         % initialize the vector of constrained dofs 
for nod = 1:size(fix_nod,2)  % for any single fixed node 
    ii = fix_nod(nod);             % ii is the number of the fixed node considered (ii=1,2,3) 
    fix_dofs = [fix_dofs 3*ii-2 3*ii-1 3*ii]; % assembly of the vector of fixed nodes, the dofs of node ii are: [3*ii-2,3*ii-1,3*ii] 
end 
free_dofs = setdiff([1:3*n],fix_dofs);  % the total number of dofs is 3*n (3 per node);  
                                                          % this function curtails the dofs that are not constrained 
% --- Loads --- %
Po = 0.01; 
f = zeros(3*n,1); 
f(3*n) = -Po; 
K_free=K(free_dofs,free_dofs); 

% --- Solution --- %
u = zeros(3*n,1);   % the function 'zeros' gives a matrix by defauls, if one of the sizes is 1, that gives a tensor 
u(free_dofs) = K(free_dofs,free_dofs)\(f(free_dofs)+fT(free_dofs)); 

% --- Plot results --- %
figure(1)
hold on 
grid on  
grid minor 
figure(2) 
hold on 
grid on 
grid minor 
for q =1:m 
    i = elem(q,1); j = elem(q,2);    
    xi = nodes(i,1); yi = nodes(i,2); zi = nodes(i,3);     % undeformed configuration  
    xj = nodes(j,1); yj = nodes(j,2); zj = nodes(j,3);  
    xid = xi+u(3*i-2); yid = yi+u(3*i-1); zid = zi+u(3*i); % deformed configuration 
    xjd = xj+u(3*j-2); yjd = yj+u(3*j-1); zjd = zj+u(3*j); 
    figure(1)  
    hold on 
    %      undeformed;               deformed;   
    plot3([xi,xj],[yi,yj],[zi,zj],'k-o',[xid,xjd],[yid,yjd],[zid,zjd],'b-o') 
    figure(2) 
    Neq = Nev(q,:); 
    edofs = [3*i-2,3*i-1,3*i,3*j-2,3*j-1,3*j]; 
    Nq = Neq*u(edofs)-kv(q)*elem(q,4)*Lv(q); 
    Nv(q) = Nq; 
    if Nq > 0  
        plot3([xi,xj],[yi,yj],[zi,zj],'r','LineWidth',Nq*400) % If the axial force is positive (tension) the color of the line is red;  
                                                                                 % The thickness of the line is proportional to the magnitude of the force; 
    elseif Nq < 0 
        plot3([xi,xj],[yi,yj],[zi,zj],'b','LineWidth',abs(Nq)*400) % If the axial force is negative (compression), the color of the line is blue; 
    end 
end 
figure(1) 
view(45,45) 
title('Structural displacements','fontsize',15) 
legend({'Udeformed cofiguration','Deformed cofiguration'},'Location','southeast') 
figure(2) 
view(45,45) 
title('Diagram of axial forces','fontsize',15); 