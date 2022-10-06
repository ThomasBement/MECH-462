%%% Cleaning memory %%% 
clear;
clc;
%%% 2D truss structure %%%

% --- Nodal coordinates --- %
%       xi yi
nodes = [0  0;          % node 1 
         1  0;          % node 2
         0 -1;          % node 3
         1 -1];         % node 4 
n = size(nodes,1);      % number of nodes 

% --- Elements --- %
%      i j EA aDT   
elem = [1 2 1 0;        % element 1
        3 4 1 0;        % element 2
        1 4 1 0;        % element 3
        3 2 1 0;        % element 4
        2 4 1 0.05];    % element 5
m = size(elem,1);       % number of elements 

% Stiffness matrix and vector of ‘thermal forces’ 
K = zeros(2*n,2*n); 
fT = zeros(2*n,1);   % initialization 

% assembly 
for q = 1:m 
    i = elem(q,1); xi = nodes(i,1); yi = nodes(i,2); % extracting the coordinates of the nodes of pertinence 
    j = elem(q,2); xj = nodes(j,1); yj = nodes(j,2); 
    L = sqrt((xj-xi)^2+(yj-yi)^2); Lv(q)=L; 
    kq = elem(q,3)/L;                % stiffness 
    kv(q) = kq; 
    qv  = ([xj;yj]-[xi;yi])/L;          % vector q 
    Ne = (kq*[-qv;qv]);             % vector Ne used to obtain the stiffness matrix and to evaluate the axial forces in the trusses 
    Nev(q,:) = Ne;                   % storing the vector Ne for element q for future use 
    Ke = 1/kq*Ne*Ne'; 
    edofs = [2*i-1,2*i,2*j-1,2*j]; 
    fT(edofs) = fT(edofs)+elem(q,4)*L*Ne;  % the vector fT relaed to element q is fT=aDTq*L*Ne, with aDTq =elem(q,4) (*) 
    K(edofs,edofs) = K(edofs,edofs)+Ke;   
end 
% Boundary conditions 
fix_dofs = [1 2 5 6];        % dofs constrained to be 0 
free_dofs = setdiff([1:2*n],fix_dofs);  % the total number of dofs is 2*n (2 per node); this function curtails the dofs that are not      
                                                          % constrained 
% Loads 
Po = 0.01; 
f = zeros(2*n,1);  
%f(4)=-Po; Add loads here
% Solution  
u = zeros(2*n,1); 
u(free_dofs) = K(free_dofs,free_dofs)\(f(free_dofs)+fT(free_dofs)); 
% Plot results 
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
    xi = nodes(i,1); yi = nodes(i,2);    % undeformed configuration  
    xj = nodes(j,1); yj = nodes(j,2); 
    xid = xi+u(2*i-1); yid = yi+u(2*i);  % deformed configuration 
    xjd = xj+u(2*j-1); yjd = yj+u(2*j); 
    figure(1)  
    %    undeformed;     deformed;   
    plot([xi,xj],[yi,yj],'k-o',[xid,xjd],[yid,yjd],'b-o') 
    figure(2) 
    Neq = Nev(q,:); 
    edofs = [2*i-1,2*i,2*j-1,2*j]; 
    Nq = Neq*u(edofs)-kv(q)*elem(q,4)*Lv(q);  % here we calculate Nq in the usual way but curtail out the artificial thermal force (*) 
    Nv(q) = Nq; 
    if Nq > 0  
        plot([xi,xj],[yi,yj],'r','LineWidth',Nq*400)    % if the axial force is positive, i.e. it creates tension, then the color is red 
                                                                        % the thickness of the line is proportional to the magnitude of the force 
    elseif Nq < 0 
        plot([xi,xj],[yi,yj],'b','LineWidth',abs(Nq)*400) % if the axial force is negative, the color of the line is blue; 
    end 
end 
figure(1) 
title('Structural displacements','fontsize',15) 
legend({'Udeformed cofiguration','Deformed cofiguration'},'Location','southeast') 
figure(2) 
title('Diagram of axial forces','fontsize',15); 
 