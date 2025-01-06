% I was not assigned parameters in project_values excel file, 
% so I used the parameters in the first row of the file

% steady-state temperature in our cylinder
%% parameters, initializing and setting boundary conditions
radius = 0.2;           % radius of the cylinder
height = 2.5;           % height of the cylinder
lambda = 406;           % conductivity of the material
T0 = -25;               % ambient temperature
T1 = 215;               % bottom end

T_initial = 100;        % initial temperature inside the cylinder
                        % this can be any value, because it does not affect
                        % the steady state response, but it can change the
                        % number of iterations to get to the steady state.
                        % Steady-state response is only affected by the
                        % boundary condition

dr = 0.01;              % by reducing this value, the resolution of the
                        % grid network increases but in return 
                        % the computation takes longer

dz = dr;                % using eaual steps for r and z directions

tolerance = 1e-4;       % tolerance between two consecutive values
                        % for temperature at a grid point. By reducing this
                        % value, the accuracy of the temperature value at
                        % each points increase but

M = radius / dr + 1;    % number of radial points
N = height / dz + 1;    % number of points in z direction
r = 0:dr:radius;
z = 0:dz:height;

% making the grid points matrix and 
% setting the initial temperature of the cylinder
T = T_initial * ones(M,N);

% setting the boundary conditions
T(:,1) = T1;            % bottom of the cylinder at T1
T(:,end) = T0;          % top of the cylinder at T0
T(end,:) = T0;          % side of the cylinder at T0

max_iteration = 10000;     % maximum number of iterations 
                           % stops the computation, in case the model 
                           % does not achieve the desired value of
                           % tolerance before 10000 iterations

%% The Gauss-Seidel method  

for n = 1:max_iteration
    T_old = T;
    for i = 2:M
        for j = 2:N
            % formula has been discussed in the report

            % ef is an auxiliary function defined to help with
            % handling boundary conditions during calculations

            % As it is obvious from the formula, we use the T matrix
            % to update the value of the temperature at different points
            % So, this is the the Gauss-Seidel method. Instead, if
            % we wanted to use the Jacobi method, we had to
            % use T_old as the argument to the each ef function evaluation
            % To use Jacobi method instead, uncomment the next section and
            % comment this one.

            T(i,j) = (ef(T,i+1,j)/dr^2 ... 
                + ef(T,i-1,j)/dr^2 ... 
                + (1/2/dr/((i-1)*dr))*(ef(T,i+1,j)-ef(T,i-1,j)) ...
                + ef(T,i,j+1)/dz^2 ... 
                + ef(T,i,j-1)/dz^2) ...
                / (2/dr^2 + 2/dz^2);
        end
    end
    % mantaining boundary conditions in case they have changed
    T(:,1) = T1;            % bottom of the cylinder at T1
    T(:,end) = T0;          % top of the cylinder at T0
    T(end,:) = T0;          % side of the cylinder at T0
    T(1,:) = T(2,:);        % otherwise temperature won't change at the z 
                            % axis
    if max(max(abs(T_old-T))) < tolerance
        fprintf('The number of iterations: %d iterations\n',n)
        break;
    end
end

%% The Jacobi method  
%{
for n = 1:max_iteration
    T_old = T;
    for i = 2:M
        for j = 2:N
            % Jacobi method
            T(i,j) = (ef(T_old,i+1,j)/dr^2 ... 
                + ef(T_old,i-1,j)/dr^2 ... 
                + (1/2/dr/((i-1)*dr))*(ef(T_old,i+1,j)-ef(T_old,i-1,j)) ...
                + ef(T_old,i,j+1)/dz^2 ... 
                + ef(T_old,i,j-1)/dz^2) ...
                / (2/dr^2 + 2/dz^2);
        end
    end
    % mantaining boundary conditions in case they have changed
    T(:,1) = T1;            % bottom of the cylinder at T1
    T(:,end) = T0;          % top of the cylinder at T0
    T(end,:) = T0;          % side of the cylinder at T0
    T(1,:) = T(2,:);        % otherwise temperature won't change at the z 
                            % axis
    if max(max(abs(T_old-T))) < tolerance
        fprintf('The number of iterations: %d iterations\n',n)
        break;
    end
end
%}

%% plotting the steady-state temperature
[R,Z] = meshgrid(r,z);
contourf(R,Z,T')
colorbar
title("Steady-state Temperature")
xlabel("r")
ylabel("z")

%% calculating tmperature gradient
rComponent = zeros(M,N);
zComponent = zeros(M,N);

for i = 1:M
    for j = 1:N
        % T now represnts the steady state temperature

        % the divergence of temperature gradient is zero everywhere
        % becasue we are in the steady state

        % due to large range of values for temperature graident 
        % at different points, this vector field has not been plotted
        % but the component values can be accessed as worksapce variables

        rComponent(i,j) = (ef(T,i+1,j) - ef(T,i-1,j))/2/dr;
        zComponent(i,j) = (ef(T,i,j+1) - ef(T,i,j-1))/2/dz;
    end
end

%% an auxiliary function to handle boundaries
% "ef" stands for entry finder
function output = ef(T,i,j)   
    [row,column] = size(T);
    if (i >= 1 && i <= row) && (j >= 1 && j <= column)
        output = T(i,j);
    end
    if i == 0
        output = T(1,j);
    end
    if j == 0
        output = T(i,1);
    end
    if i == row + 1
        output = T(row,j);
    end
    if j == column + 1
        output = T(i,column);
    end
end
    

