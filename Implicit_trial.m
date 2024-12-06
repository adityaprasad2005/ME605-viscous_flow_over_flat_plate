clc
clear all
close all

% Parameters
L = 1;
H = 0.1;
Re = 1e4;
U_inf = 1;
nu = U_inf * L / Re;
tol = 1e-6;
max_iter = 10000;
nx = 1000;
ny = 1000;
dx = L / (nx - 1);
dy = H / (ny - 1);

% Initialize u and v velocity fields
u = U_inf.*ones(ny, nx);
v = zeros(ny, nx);

% Boundary Conditions
% Left boundary
u(:, 1) = U_inf;
% Bottom boundary
u(1,:) = 0;
% Top boundary
u(end, :) = U_inf;


for j = 2:nx
    u(:,j) = u(:,j-1);
    no_convergence = true;
    iter = 0;
    while (no_convergence)
        u_old = u(:,j);
        for i = 2:ny-1
            u(i,j) = solve_u(i, j, u, v, nu, dx, dy);
        end
        residual = sqrt(sum((u(:,j) - u_old).^2));
        iter = iter + 1;
        if(residual<tol || iter>=max_iter)
            no_convergence = false;
        end
    end

    v(:,j) = v(:,j-1);
    no_convergence = true;
    iter = 0;
    while (no_convergence)
        v_old = v(:,j);
        for i = 2:ny-1
            v(i,j) = solve_v(i, j, u, v, dx, dy);
        end
        residual = sqrt(sum((v(:,j) - v_old).^2));
        iter = iter + 1;
        if(residual<tol || iter>=max_iter)
            no_convergence = false;
        end
    end
end

lis = [];
lis = [lis,0];
lis2 = [];
lis2 = [lis2,1];

for j = 2:nx
    for i = 2:ny-1
        if u(i,j)>= 0.9999*U_inf
            lis = [lis,i];
            lis2 = [lis2,j];
            break;
        end
    end
end

x_arr = linspace(0, L, nx);
y_arr = linspace(0, H, ny);
[X,Y] = meshgrid(x_arr,y_arr);

% Plotting u-velocity contour
figure;
contourf(X,Y, u, 20);
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
contourf(linspace(0, L, nx), linspace(0, H, ny), v);
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

disp('ended');

% Solver for u-velocity
function f = solve_u(i, j, u, v, nu, dx, dy)
    e = u(i,j-1);
    c = u(i+1,j);
    ai = u(i-1,j);
    f = ((e^2/dx) - c*(v(i,j-1)/(2*dy)-nu/(dy^2)) + ai*(v(i,j-1)/(2*dy)+nu/(dy^2)))/((e/dx)+(2*nu/(dy^2)));
    % if f<0
    %     f=0;
    % end
end

% Solver for v-velocity
function f = solve_v(i, j, u, v, dx, dy)
    ai = v(i-1,j);
    f = ai + (u(i,j)-u(i,j-1))*dy/dx;
    % if f<0
    %     f=0;
    % end
end