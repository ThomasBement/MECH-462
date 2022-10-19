%%% Dynamic analysis of a 3D structure  %%%
%%% Compares implicit versus explicit   %%%

clc;
clear;  

% --- NODES --- % 
%       [xi yi zi] -> assume: element length a = 1
nodes = [0 0 0;                             % node 1 
        sqrt(3)/2 0.5 0;                    % node 2    
        0 1 0;                              % node 3 
        0.5/sqrt(3) 0.5 sqrt(2/3)];         % node 4 
n = size(nodes,1);                          % number of nodes 

% --- ELEMENTS --- %
%      [i j EA] -> assume: nominal EAo = 1 (**) 
elem = [1 4 1;                              % element 1 
        2 4 1;                              % element 2 
        3 4 1];                             % element 3 
m = size(elem,1);                           % number of elements 

% --- STIFFNESS MATRIX --- %
K = zeros(3*n,3*n);  % initialization;

% --- ASSEMBLY --- %
for q = 1:m 
    i = elem(q,1); xi = nodes(i,1); yi = nodes(i,2); zi = nodes(i,3); % extracting the coordinates of the nodes of pertinence 
    j = elem(q,2); xj = nodes(j,1); yj = nodes(j,2); zj = nodes(j,3); 
    L = norm([xj;yj;zj]-[xi;yi;zi]);  Lv(q)=L; 
    kq = elem(q,3)/L;  % stiffness 
    qv = ([xj;yj;zj]-[xi;yi;zi])/L; 
    Ne = (kq*[-qv;qv]); 
    Nev(q,:) = Ne; 
    Ke = 1/kq*Ne*Ne'; 
    edofs = [3*i-2,3*i-1,3*i,3*j-2,3*j-1,3*j]; 
    K(edofs,edofs) = K(edofs,edofs)+Ke;  
end 

% --- BOUNDARY CONDITIONS --- %
fix_nod = [1 2 3];  % nodes that are fixed 
fix_dofs = [];         % initialize the vector of constrained dofs 
for nod = 1:size(fix_nod,2)  % for any single fixed node 
    ii = fix_nod(nod);             % ii is the number of the fixed node considered (ii=1,2,3) 
    fix_dofs = [fix_dofs 3*ii-2 3*ii-1 3*ii]; % assembly of the vector of fixed nodes, the dofs of node ii are: [3*ii-2,3*ii-1,3*ii] 
end 
free_dofs = setdiff([1:3*n],fix_dofs);  % the total number of dofs is 3*n (3 per node);  
                                                        % this function curtails the dofs that are not constrained 
nf = size(free_dofs,2);  % number of free dofs 

% --- INITIAL CONDITIONS --- %
u0 = [0.1;0;0];     % initial displacements of node 4 (the only one free) 
du0 = [0;0;0];      % initial velocity of node 4 (the only one free) 

% --- MASS --- % 
M = eye(3*n);       % eye creates an identity matrix, 
                    % we take the mass m = 1 (**) 

% --- IMPLICIT ANALYSIS --- %
B = K(free_dofs,free_dofs)/M(free_dofs,free_dofs); 
[V,Bdiag] = eig(B);                 % here V is the matrix including 
                                    % the eigenvectors, hence R', 
                                    % Bdiag is the diagonal matrix of B 
omegas = sqrt(Bdiag*ones(nf,1));    % this is a vector which components 
                                    % are the natural frequencies omega 

% The characteristic time t* scales with the period T=2*pi/omega  
tf = 20;                    % maximum simulation time 
nt_steps = 100;             % number of time steps for visualization; 
                            % For the implicit analysis this does 
                            % not affect the accuracy 
Dt = tf/nt_steps;           % time step for visualization 
tv_imp = [0:Dt:tf]; 
Abv = V'*u0;                % this is the coefficient of cos(omega*t) V'=R 
Av = sqrt(Bdiag)\V'*du0;    % this is the coefficient of sin(omega*t)
for i = 1:size(tv_imp,2) 
    % Store position vector
    uv_imp(i,:) = V*(Av.*sin(omegas*tv_imp(i))+Abv.*cos(omegas*tv_imp(i)));
    % Store velocity vector
    duv_imp(i,:) = V*(Av.*omegas.*cos(omegas*tv_imp(i))-Abv.*omegas.*sin(omegas*tv_imp(i)));
    % Store acceleration vector
    dduv_imp(i,:) = -V*(Av.*omegas.*omegas.*sin(omegas*tv_imp(i))+Abv.*omegas.*omegas.*cos(omegas*tv_imp(i)));
end 

figure(1) 
hold on 
plot(tv_imp,uv_imp(:,1),"o",'DisplayName','Implicit: Position');
plot(tv_imp,duv_imp(:,1),"o",'DisplayName','Implicit: Velocity');
plot(tv_imp,dduv_imp(:,1),"o",'DisplayName','Implicit: Acceleration');

% --- EXPLICIT ANALYSIS --- %
t = 0; 
omax = max(omegas); 
dtmax = 1/(2*omax);     % This is the maximum time step we discussed in class; it works best when the initial displacement is 0 ... 
dt = dtmax;      
u = u0;                     % Initializing position
du = du0;                   % Initializing velocity
ddu = -B*u0*dt;             % Initializing acceleration
tv_exp = [0:dt:tf];         % Vector of time-steps for explicit analysis 
for i = 1:size(tv_exp,2) 
    uv_exp(i,:) = u;        % Storing position vector
    duv_exp(i,:) = du;      % Storing velocity vector
    dduv_exp(i,:) = ddu;    % Storing acceleration vector
    
    u = u+du*dt;            % Iterating position
    ddu = -B*u;             % Iterating acceleration
    du = du+ddu*dt;         % Iterating velocity
    t = t+dt;               % Iterating time
end 
figure(1) 
plot(tv_exp,uv_exp(:,1),'DisplayName','Explicit: Position')
plot(tv_exp,duv_exp(:,1),'DisplayName','Explicit: Velocity')
plot(tv_exp,dduv_exp(:,1),'DisplayName','Explicit: Acceleration')
legend