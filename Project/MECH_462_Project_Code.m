%%% MECH-462 Project Code %%% 
clear;  
%-Problem Definition--------%------------------------------------
a       = 1.0;              % Punch Radius                      [m]
delta   = -1e-5;            % Punch Displacement                [m]
R       = 20*a;             % Sample Radius                     [m]
H       = R;                % Sample Height                     [m]
E_s     = 70e9;             % Elastic Modulus                   [Pa]
V_s     = 0.3;              % Poissons Ratio                    [N.a.]
%-Plot parameters-----------%------------------------------------
nms     = 2;                % Node Marker Size
scale   = 1;              % Displacements Scale Factor
%-Material Definition-------%------------------------------------
E       = E_s;              % Elastic Modulus 
v       = V_s;              % Poisson's Ratio 
%-Structure Definition------%------------------------------------
Lx      = R;                % Length in X 
Ly      = H;                % Length in Y
%-Mesh----------------------%------------------------------------
nx      = 50;               % Number of X Elements
ny      = 50;               % Number of Y Elements
m       = nx*ny;            % Total Number of Elements
Dx      = Lx/nx;            % Size of the in X 
Dy      = Ly/ny;            % Size of the in Y
%-Nodes---------------------%------------------------------------
N       = (nx+1)*(ny+1);    % Total Number of Nodes
ix      = [1:nx+1];         % X Component of Node Grid
iy      = [1:ny+1];         % Y Component of Node Grid
x       = (ix-1)*Dx;        % X Values of Node Grid
y       = (iy-1)*Dy;        % Y Values of Node Grid
[xn,yn] = meshgrid(x,y);    % Node Grid
%-Stiffness Matrix----------%------------------------------------
K   =   zeros(2*N,2*N);     % Initializing K Matrix
Em  =   [(1-v) v    0;      % 
           v (1-v)  0;      %
           0   0 (1-2*v)/2];% 
Em  = Em*E/((1+v)*(1-2*v)); % Plane Strain Elasticity Matrix
T   =   [1,0,0,0;           % Derivative Mapping Matrix
         0,0,0,1;           %
         0,1,1,0];          %  
g(1,:) = 1/sqrt(3)*[-1,-1]; % Gauss Points Matrix
g(2,:) = 1/sqrt(3)*[1,-1];  %
g(3,:) = 1/sqrt(3)*[-1,1];  %
g(4,:) = 1/sqrt(3)*[1,1];   %
%-Assembly------------------%-----------------------------------%------------------------------------
for ex = 1:nx                                                   %
    for ey = 1:ny                                               %
        i = (ey-1)*(nx+1)+ex;                                   % Index for bottom row of element
        j = i+nx+1;                                             % Index for top row of element
        x1 = xn(ey,ex); x2 = xn(ey,ex+1); x3 = x1; x4 = x2;     % X coordinates for elemnt nodes
        y1 = yn(ey,ex); y2 = y1; y3 = yn(ey+1,ex); y4 = y3;     % Y coordinates for elemnt nodes
        Kq     = 0;                                             % Initializing Kq element
        for gg = 1:4                                            % Integrate using Gauss points 
            xx = g(gg,1);                                       %
            yy = g(gg,2);                                       %
            Nxx = [-0.25*(1-yy), 0.25*(1-yy), -0.25*(1+yy), 0.25*(1+yy)];                       % Shape function 
            Nyy = [-0.25*(1-xx), -0.25*(1+xx), 0.25*(1-xx), 0.25*(1+xx)];                       % Shape function
            J = [Nxx*[x1;x2;x3;x4], Nxx*[y1;y2;y3;y4]; Nyy*[x1;x2;x3;x4], Nyy*[y1;y2;y3;y4]];   % Jacobian
            dJ = det(J);                                        % Determinate of Jacobian
            iJ = inv(J);                                        % Inverse of Jacobian
            A = T*[iJ,zeros(2,2);zeros(2,2),iJ];                % A matrix
            G = [Nxx(1) 0  Nxx(2) 0  Nxx(3)  0   Nxx(4) 0;      % G matrix
                 Nyy(1) 0  Nyy(2) 0  Nyy(3)  0   Nyy(4) 0;      %
                   0  Nxx(1) 0  Nxx(2) 0   Nxx(3)  0  Nxx(4);   %
                   0  Nyy(1) 0  Nyy(2) 0   Nyy(3)  0  Nyy(4)];  %
            B = A*G;                                            % B matrix
            Kq = Kq + B'*Em*B*dJ;                               % Iterating Kq element
        end                                                     %
        edofs = [2*i-1,2*i,2*(i+1)-1,2*(i+1),2*j-1,2*j,2*(j+1)-1,2*(j+1)];  % DOFs for Kq with respect to K matrix
        K(edofs,edofs) = K(edofs,edofs)+Kq;                                 % Fill in K matrix with Kq
    end                                                         %
end                                                             %
%-Boundary conditions-----------%-------------------------------%----
u = zeros(2*N,1);               % Intializing u vector
fix_dofs = [];                  % Intializing fixed DOFs vector 
%-Loads-------------------------%------------------------------------ 
f = zeros(2*N,1);               % Initializing force vector
%-Prescribed Displacements------%---------------%--------------------
for i = 1:N                                     %
    i_x = 2*i-1;                                %
    i_y = 2*i;                                  %
    row = floor((i-1)/(nx+1));                  %
    col = mod(i-1, (ny+1));                     %
    %-Fixed Bottom Row--------------------------%------------------------------------
    if (row == 0/Dy)                            %
        u(i_x) = 0;                             %
        u(i_y) = 0;                             %
        fix_dofs = [fix_dofs, i_x];             %
        fix_dofs = [fix_dofs, i_y];             %
    %-Symmetry Condition------------------------%------------------------------------
    elseif (col == 0/Dy) && ~(row == H/Dy)      %
        u(i_x) = 0;                             %
        fix_dofs = [fix_dofs, i_x];             %
    %-Perscribed Displacements for Punch--------%------------------------------------
    elseif (row == H/Dy) && (col <= a/Dx)       %
        u(i_x) = 0;                             %
        u(i_y) = delta;                         %
        fix_dofs = [fix_dofs, i_x];             %
        fix_dofs = [fix_dofs, i_y];             %
    end                                         %
end                                             %
iypr = fix_dofs(2:2:end);                       % PD Y 
ixpr = fix_dofs(1:2:end);                       % PD X
%-Solution--------------------------------------%------------------------------------
free_dofs = setdiff([1:2*N],fix_dofs);          %
freeNP_dofs = setdiff(free_dofs,fix_dofs);      % free dofs with Non Prescribed displacements 
K_ff = K(freeNP_dofs,freeNP_dofs);              %
K_fp = K(freeNP_dofs,fix_dofs);                 %
u(freeNP_dofs) = K_ff\(f(freeNP_dofs)-K_fp*u(fix_dofs));    % Solve problem
writematrix(u,'matlab_u_vect.txt');             % Write displacement vector
%-Plot results----------------------------------%------------------------------------
figure(1) 
hold on; grid on; grid minor 
title('Nodes','fontsize',15) 
%plot(xn(iypr,ixpr),yn(iypr,ixpr),'bs','Markersize',6) 
plot(xn,yn,'ko','Markersize',3) 
legend('prescribed displacements nodes','all nodes','Location','NorthEastOutside') 
 
axis equal; 
%%% 
figure(2) 
title('Displacements','fontsize',15) 
hold on; grid on; grid minor 
grey = 0.6*[1 1 1]; 
for ey = 1:ny 
    for ex = 1:nx 
        % Position information
        x1 = xn(ey,ex); x2 = xn(ey,ex+1);   x3 = x1;            x4 = x2; 
        y1 = yn(ey,ex); y2 = y1;            y3 = yn(ey+1,ex);   y4 = y3;  
        % Indexing
        i = (ey-1)*(nx+1)+ex; 
        j = i+nx+1;
        % Displacement Information
        u1 = scale*u(2*i-1,1);      v1 = scale*u(2*i,1);   
        u2 = scale*u(2*(i+1)-1,1);  v2 = scale*u(2*(i+1),1);  
        u3 = scale*u(2*j-1,1);      v3 = scale*u(2*j,1);  
        u4 = scale*u(2*(j+1)-1,1);  v4 = scale*u(2*(j+1),1); 
        % Plotting deformation grid
        plot([x1,x2],[y1,y2],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey) 
        plot([x1,x3],[y1,y3],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey) 
        plot([x2,x4],[y2,y4],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey) 
        plot([x3,x4],[y3,y4],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey) 
        plot([x1+u1,x2+u2],[y1+v1,y2+v2],'b-o','MarkerSize',nms,'markerfacecolor','b') 
        plot([x1+u1,x3+u3],[y1+v1,y3+v3],'b-o','MarkerSize',nms,'markerfacecolor','b') 
        plot([x2+u2,x4+u4],[y2+v2,y4+v4],'b-o','MarkerSize',nms,'markerfacecolor','b') 
        plot([x3+u3,x4+u4],[y3+v3,y4+v4],'b-o','MarkerSize',nms,'markerfacecolor','b') 
        % Strain information
        ee11 = 0.5*((u2-u1)/(x2-x1)+(u4-u3)/(x4-x3)); 
        ee22 = 0.5*((v3-v1)/(y3-y1)+(v4-v2)/(y4-y2)); 
        ee12 = 0.5*((u3-u1)/(y3-y1)+(u4-u2)/(y4-y2)+(v2-v1)/(x2-x1)+(v4-v3)/(x4-x3)); 
        eev(ey,ex,:) = [ee11;ee22;ee12];
        % Stress information
        sv(ey,ex,:) = Em*[ee11;ee22;ee12];
        % Strain energy
        Uv(ey,ex) = 0.5*[ee11,ee22,ee12]*Em*[ee11;ee22;ee12];
        % Position information
        xc(ey,ex) = 0.25*(x1+x2+x3+x4); 
        yc(ey,ex) = 0.25*(y1+y2+y3+y4); 
    end 
end

axis equal
figure(40)
plot(xc(end,:), sv(end,:,2))



axis equal 
for splt = 1:3 
    figure(10+splt)
    title(sprintf('Strain %i', splt),'fontsize',15) 
    contourf(xc,yc,sv(:,:,splt)) 
    colorbar 
    axis equal 
end 
for eplt = 1:3
    figure(20+eplt)
    title(sprintf('Stress %i', eplt),'fontsize',15) 
    contourf(xc,yc,eev(:,:,eplt)) 
    colorbar 
    axis equal 
end 
figure(30)
title('Strain Energy','fontsize',15) 
contourf(xc,yc,Uv) 
colorbar 
axis equal 