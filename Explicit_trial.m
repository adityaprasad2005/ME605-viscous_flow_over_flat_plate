clc
clear all
close all

% Parameters
L = 1;
H = 0.2;
Re = 1e4;
U_inf = 1;
nu = U_inf * L / Re;

nx = 10000;
ny = 200;
dx = L / (nx - 1);
dy = H / (ny - 1);

% Initialize u and v velocity fields
u = zeros(ny, nx);
v = zeros(ny, nx);

% Boundary Conditions
% Left boundary
u(:, 1) = U_inf;
% Top boundary
u(end, :) = U_inf;

lis = [];
lis = [lis,0];
lis2 = [];
lis2 = [lis2,1];

for j = 2:nx
    flag = false;
    for i = 2:ny-1
        u(i, j) = solve_u(i, j, u, v, nu, dx, dy);
        if flag==false && u(i,j)>=0.99*U_inf
            flag = true;
            lis = [lis,i];
            lis2 = [lis2,j];
        end
        v(i, j) = solve_v(i, j, u, v, dx, dy);
    end
end

x_arr = linspace(0, L, nx);
y_arr = linspace(0, H, ny);


% Create the matrix with ny rows, each row being x_arr
x = repmat(x_arr, ny, 1);

y = repmat(y_arr', 1, nx);  % Transpose y_arr to make it a column vector

[X,Y] = meshgrid(x_arr,y_arr);

% Plotting u-velocity contour
figure;
contourf(X,Y, u, 10);
colorbar;
new_lis = (H/ny) .* lis;
hold on;
new_lis2 = (L/nx) .* lis2;
plot(new_lis2, new_lis);
title('Contour of u-velocity');
xlabel('x');
ylabel('y');

% Plotting v-velocity contour
figure;
contourf(X,Y,v,50);
colorbar;
title('Contour of v-velocity');
xlabel('x');
ylabel('y');

% Compute and plot normalized x-velocity (F') as a function of similarity variable (eta)
eta = linspace(0, H * sqrt(U_inf / nu / L), ny);
F_prime = u(:, end) / U_inf;

figure;
plot(eta, F_prime, '-o');
title('Normalized x-velocity F''(\eta) vs Similarity Variable \eta');
xlabel('Similarity Variable \eta');
ylabel('Normalized x-velocity F''(\eta)');

% Compute and plot boundary layer thickness (delta)
delta_x = 4.91 .* sqrt(nu*linspace(0, L, nx)/U_inf);
figure;
plot(linspace(0, L, nx), delta_x);
new_lis = (H/ny) .* lis;
hold on;
new_lis2 = (L/nx) .* lis2;
plot(new_lis2, new_lis);
title('Boundary Layer Thickness \delta vs x');
xlabel('x');
ylabel('Boundary Layer Thickness \delta');


% Solver for u-velocity
function e = solve_u(i, j, u, v, nu, dx, dy)
    a = u(i+1, j-1);
    d = u(i, j-1);
    g = u(i-1, j-1);
    e = d + nu * (a - 2*d + g) / d * dx / (dy^2) - (a - g) / 2 * v(i, j-1) / d * dx / dy;
    % if d==0 || e<0
    %     e=0;
    % end
end

% Solver for v-velocity
function e = solve_v(i, j, u, v, dx, dy)
    h = v(i-1, j);
    e = h - 0.5 * (u(i,j) - u(i,j-1) + u(i-1,j) - u(i-1,j-1)) * dy / dx;
    % if e<0
    %     e=0;
    % end
end

