
close, clear, clc

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp = strjoin(splittedstr(1:end-1), '/'); 
addpath(genpath(sprintf('%s/Functions', pp)))

fig_folder = sprintf('%s/Figures', current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end

%% Pecura & Caroll Configuration

% Lorenz System driven by x, y and z variables respectively.

% Pecora-Caroll Configuration - Lorenz drive and response with conditional
% Lyapunov of the system shown in the paper 'The synchronization of chaotic
% systems' by Boccaletti et.al. 2002.

%% Drive System >>> Lorenz
%  Response >> Lorenz which has been driven by x-variable of drive.
%  Replika >> Lorenz which has been driven by x-variable of drive.

% Conditional Lyapunov Exp >> Lyapunov exponents for the subsystem of
% response, its negative value indicates the sync of the orbits of response
% and its replika. But not guarantee that the they will stay in sync
% forever.

x0 = 100*rand(1,7); % 7 initial conditions, 3 drive systems components, 2 for the replika of response y & z, 2 for response.
t0 = 0;
tf = 1e3;
k = 1e4;

len_of_deviation_vector = 4; % Response system 2d, Jacobian matrix has 2-by-2 entries. 

[tvalues_drivenby_x, conditional_LYA_drivenby_x, Orbits_drivenby_x] = LyapunovforODEs_defined_length_deviation(x0,t0, tf, k, @x_driver_response_Lorenz, len_of_deviation_vector);
% Orbits have 7 columns, 1:3 >> drive parameters, 4:5 >> Response, 6:7 >>
% Replika of response % Response variables driven by x of driver.

[tvalues_drivenby_y, conditional_LYA_drivenby_y, Orbits_drivenby_y] = LyapunovforODEs_defined_length_deviation(x0,t0, tf, k, @y_driver_response_Lorenz, len_of_deviation_vector);
% Orbits have 7 columns, 1:3 >> drive parameters, 4:5 >> Response, 6:7 >>
% Replika of response % Response variables driven by x of driver.

%
x0 = ones(1,7);
[tvalues_drivenby_z, conditional_LYA_drivenby_z, Orbits_drivenby_z] = LyapunovforODEs_defined_length_deviation(x0,t0, tf, k, @z_driver_response_Lorenz, len_of_deviation_vector);
% Orbits have 7 columns, 1:3 >> drive parameters, 4:5 >> Response, 6:7 >>
% Replika of response % Response variables driven by x of driver.
%
% Absolute error between y,z and y', z' to show that it loses sync.
e1 = abs(Orbits_drivenby_z(:,4)-Orbits_drivenby_z(:,6)); % Absolute error of x and x'.
e2 = abs(Orbits_drivenby_z(:,5)-Orbits_drivenby_z(:,7)); % Absolute error of y and y'.
%%

figure
% subplot(2,2,1)
% plot3(Orbits_drivenby_x(:,1), Orbits_drivenby_x(:,2), Orbits_drivenby_x(:,3), '.')
% title('Drive System', 'Interpreter','latex', 'FontSize',20)

subplot(2,2,1)
plot(tvalues_drivenby_x(1:5e1), Orbits_drivenby_x(1:5e1,4), 'LineWidth',2, 'LineStyle',':') % response y
hold on
plot(tvalues_drivenby_x(1:5e1), Orbits_drivenby_x(1:5e1,6), 'LineWidth',2)  % replika of y.
ax = gca;
ax.FontSize = 30; 
title('(a)', 'FontSize',30, 'Interpreter','latex')
% ylabel("$y_{1}, y'_{1}$", 'FontSize',30, 'Interpreter','latex')
xlabel('t', 'FontSize',30, 'Interpreter','latex')
ylim([-100, 50])
yticks([-100, -50, 0, 50])
legend('$y_{1}$',"$y'_{1}$", 'FontSize', 30, 'Interpreter','latex','Location', 'best')


subplot(2,2,3)
plot(tvalues_drivenby_x(1:5e1), Orbits_drivenby_x(1:5e1,5), 'LineWidth',2, 'LineStyle',':') % response z
hold on
plot(tvalues_drivenby_x(1:5e1), Orbits_drivenby_x(1:5e1,7), 'LineWidth',2)  % replika of z.
ax = gca;
ax.FontSize = 30; 
title('(b)', 'FontSize',30, 'Interpreter','latex')
xlabel('t', 'FontSize',30, 'Interpreter','latex')
% ylabel("$z_{1}, z'_{1}$", 'FontSize',30, 'Interpreter','latex')
ylim([0, 100])
yticks([0, 50, 100])
legend('$z_{1}$',"$z'_{1}$", 'FontSize', 30, 'Interpreter','latex','Location', 'best')

subplot(2,2,[2,4])
plot(tvalues_drivenby_z(1:1e1),e1(1:1e1), 'o--', 'LineWidth',2)
hold on
plot(tvalues_drivenby_z(1:1e1), e2(1:1e1),  '*--', 'LineWidth',2)
ax = gca;
ax.FontSize = 30; 
title('(c)', 'FontSize',30, 'Interpreter','latex')
xlabel('t', 'FontSize',30, 'Interpreter','latex')
% ylabel("$x_{1}, x'_{1} \\ y_{1}, y'_{1}$", 'FontSize',30, 'Interpreter','latex')
legend("$|x_{1} - x'_{1}|$", "$|y_{1} - y'_{1}|$", 'FontSize', 30, 'Interpreter','latex','Location', 'best')
%%

saveas(gca, sprintf('%s/PC_completesync.fig', fig_folder))
saveas(gca, sprintf('%s/PC_completesync', fig_folder), 'epsc')
close

%% table for conditional MLE

Driver = {'x'; 'y'; 'z'};
MLE1 = [conditional_LYA_drivenby_x(1); conditional_LYA_drivenby_y(1); conditional_LYA_drivenby_z(1)];
MLE2 = [conditional_LYA_drivenby_x(2); conditional_LYA_drivenby_y(2); conditional_LYA_drivenby_z(2)];

table_cond_MLE = table(Driver, MLE1, MLE2);
writetable(table_cond_MLE, sprintf('%s/ConditionalMLE.txt', fig_folder))

%% Bidirectionally Coupling Lorenz System  : All coupled
% Maximum Lyapunov Exponent for uncoupled system is 1.37. The oscillators
% becoming synchronized if the coupling strength is higher than 1.37/2
% (Boccaletti et. al. 2002), all variables coupled.

%   clear, close, clc
% 
x0 = 10*rand(1,6);
t0 = 0;
tf = 1e3; 
k = 1e5;  
coupling_strength = 2; 

sigma = 10; 
rho = 28;  
beta = 8/3;

Adj = [0, 1; 1 0];


[tvalues4, Lyapunov_K_2, Orbits4] = Lyapunovfor_all_coupledLorenz(x0, t0, tf, k, coupling_strength, sigma, rho, beta, Adj);


coupling_strength = 0; 
[tvalues, Lyapunov_spectrum, ~] = Lyapunovfor_all_coupledLorenz(x0, t0, tf, k, coupling_strength, sigma, rho, beta, Adj);

  %

%% All coupled Lorenz system data generation for 2 oscillators

        y0 = 10*rand(1,6);
        Adj = [0 1; 1 0];
        n=50;
        coupling_strength = linspace(0,2,n);

        sigma = 10; 
        rho = 28;  
        beta = 8/3;
        
        t0=0;
        tf=1e3;
        k=1e5;
    
         mean_error = zeros(n,1);

   for i=1:n
        tspan= linspace(t0,tf,k+1);
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) all_coupled_Lorenz(t,y, Adj, coupling_strength(i), sigma, rho, beta), tspan, y0);
        mean_error(i) = mean(sqrt((Y(:,1) - Y(:,2)).^2 + (Y(:,3) - Y(:,4)).^2 + (Y(:,5) - Y(:,6)).^2));
        sprintf('%d out of %d', i, n)
   end

%% Sync of the variables of Lorenz
% Coupling streng is 2.

figure
subplot(3,2,1)
plot(T(1:3e2),Y(1:3e2,1), 'LineStyle',':', 'LineWidth',2)
hold on
plot(T(1:3e2),Y(1:3e2,2), 'LineWidth',2)
% xlabel('T', 'Interpreter','latex','FontSize', 20)
% title("Sync of x variable", 'Interpreter','latex','FontSize', 20)
legend('$x_{1}$',"$x_{2}$", 'interpreter', 'latex', 'FontSize',30)
ax=gca;
ax.FontSize = 30; 
ax.XTickLabel = '';
title('(a)', 'FontSize',30, 'Interpreter','latex')

subplot(3,2,3)
plot(T(1:3e2),Y(1:3e2,3), 'LineStyle',':', 'LineWidth',2)
hold on
plot(T(1:3e2),Y(1:3e2,4), 'LineWidth',2)
% xlabel('T', 'Interpreter','latex','FontSize',20)
% title("Sync of y variable", 'Interpreter','latex','FontSize',20)
legend('$y_{1}$',"$y_{2}$", 'interpreter', 'latex', 'FontSize',30)
ax=gca;
ax.FontSize = 30; 
ax.XTickLabel = '';
title('(c)', 'FontSize',30, 'Interpreter','latex')

subplot(3,2,5)
plot(T(1:3e2),Y(1:3e2,5), 'LineStyle',':', 'LineWidth',2)
hold on
plot(T(1:3e2),Y(1:3e2,6), 'LineWidth',2)
xlabel('T', 'Interpreter','latex','FontSize',30)
% title("Sync of z variable", 'Interpreter','latex','FontSize',20)
legend('$z_{1}$',"$z_{2}$", 'interpreter', 'latex', 'FontSize',30)
ax=gca;
ax.FontSize = 30; 
title('(e)', 'FontSize',30, 'Interpreter','latex')

subplot(3,2,[2,4])
plot3(Orbits4(1e3:end-5e4,1), Orbits4(1e3:end-5e4,3), Orbits4(1e3:end-5e4,5))
xlabel('x', 'Interpreter','latex'); ylabel('y', 'Interpreter','latex');zlabel('z', 'Interpreter','latex');
ax=gca;
ax.FontSize = 30;
title('(b)', 'FontSize',30, 'Interpreter','latex')

subplot(3,2,[6])
plot(coupling_strength/Lyapunov_spectrum(1), mean_error, 'ko--', 'LineWidth',2)
xlabel('$\frac{K}{MLE}$', 'Interpreter','latex', 'FontSize',30)
ylabel(' $e$ ', 'Interpreter','latex', 'FontSize',30)
% title('Error function of $\frac{K}{MLE}$', 'Interpreter','latex', 'FontSize',20)
ax=gca;
ax.FontSize = 30;
title('(d)', 'FontSize',30, 'Interpreter','latex')

%%
ax.YTick = [0, 4, 8, 12,16];
ax.YTickLabel = {'0', '4', '8', '12', '16'}; 
%%
saveas(gca, sprintf('%s/sync_of_Lorenz.fig', fig_folder))
saveas(gca, sprintf('%s/sync_of_Lorenz', fig_folder), 'epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logistic Map - CS with 5 nodes %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data generation

A = [ 0 1 0 0 0;
      1 0 1 1 0;
      0 1 0 0 0;
      0 1 0 0 1;
      0 0 0 1 0];
% 
M = size(A,2);
x0 = rand(M,1);
niter = 1e5;
time_window = 1e2;
std_threshold = 1e-3;
max_discard = 1e4;

M = size(A,2);
x0 = rand(M,1);
alpha = 0.2;
niter = 1e5;

% [Orbits, Lya_overtime] = Log_generator(A, alpha, x0, 4, niter, time_window, std_threshold, max_discard,1); % Orbits generator for 5-nodes.  
% writematrix(Orbits, 'Logistic_Info/Orbits.txt');
% Read data from the file
% Lya_overtime = readmatrix('Logistic_Info/MLE_0.20.txt');
%%
Orbits = readmatrix('Logistic_Info/Orbits.txt');
Lya_overtime = readmatrix('Logistic_Info/MLE_0.20.txt');
info = readmatrix('Logistic_Info/info_0.20.txt');
discarded = info(1,2);
niter = 1e5;
A = [ 0 1 0 0 0;
      1 0 1 1 0;
      0 1 0 0 0;
      0 1 0 0 1;
      0 0 0 1 0];
M = size(A,2);
%
figure

subplot(2,2,2)
plot(discarded+1:niter+discarded, cumsum(abs(Orbits(:,1)-Orbits(:,3))), 'LineWidth',2)
hold on
plot(discarded+1:niter+discarded, cumsum(abs(Orbits(:,1)-Orbits(:,4))), 'LineWidth',2)
plot(discarded+1:niter+discarded, cumsum(abs(Orbits(:,1)-Orbits(:,5))), 'LineWidth',2)
plot(discarded+1:niter+discarded, cumsum(abs(Orbits(:,2)-Orbits(:,4))), 'LineWidth',2)
plot(discarded+1:niter+discarded, cumsum(abs(Orbits(:,2)-Orbits(:,5))), 'LineWidth',2)
title('(c)', 'FontSize',30, 'Interpreter','latex')
legend('1-3', '1-4', '1-5', '2-4', '2-5', 'interpreter', 'latex', 'Location','northwest')
xlim([0,1e5+discarded])
ax=gca;
ax.FontSize = 30;
xticks([0, 5e4, 1e5])

subplot(2,2,1)
p1 = plot(graph(A), 'MarkerSize', 20, 'NodeColor','r', 'EdgeColor','k', 'LineWidth', 2, 'NodeFontSize',30);
ax = gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
axis equal
pbaspect([1 1 1])
title('(a)', 'FontSize',30, 'Interpreter','latex')

subplot(2,2,4)
semilogy(discarded+1:niter+discarded, abs(Orbits(:,1)-Orbits(:,3)), 'LineWidth',2)
% title('Absolute difference between node 1 and node 3', 'interpreter', 'Latex', 'FontSize',20), 
xlabel('Time-series length', 'interpreter', 'Latex', 'FontSize',30)
xlim([0,1e5+discarded])
% ylim([0,10])
%pbaspect([ 1 1 1])
ax=gca;
ax.FontSize = 30;
yticks([1e-15,1]);
xticks([0, 5e4, 1e5])
title('(d)', 'FontSize',30, 'Interpreter','latex')

subplot(2,2,3)
loglog(1:niter+discarded, Lya_overtime(:,1), 'LineWidth',2)
xline(discarded, 'r--', 'LineWidth',2)
% xlabel('Iteration', 'interpreter', 'Latex', 'FontSize',12)
% ylabel('MLE','interpreter', 'Latex', 'FontSize',20)
% title('Evolution of MLE', 'interpreter', 'Latex', 'FontSize',20)
xlim([0,1e5+discarded])
ylim([0,1])
ax=gca;
ax.FontSize = 30;
title('(b)', 'FontSize',30, 'Interpreter','latex')
xlabel('Time-series length', 'interpreter', 'Latex', 'FontSize',30)
ylabel('MLE', 'interpreter', 'Latex', 'FontSize',30)

%%
saveas(gca, sprintf('%s/Logistic_CS.fig', fig_folder))
saveas(gca, sprintf('%s/Logistic_CS', fig_folder), 'epsc')













%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = x_driver_response_Lorenz(t,x0)
% Parameters
sigma = 16; b = 4; r = 45.92;

% Driver system (Lorenz)
x = x0(1); y = x0(2); z = x0(3);

x_dot = sigma*(y-x);
y_dot = -x*z+r*x-y;
z_dot = x*y-b*z;

driver = [x_dot; y_dot; z_dot];

% Replika of Response
y2 = x0(6); z2 = x0(7);

y_dot2 = -x*z2 + r*x-y2;
z_dot2 = x*y2-b*z2;

% Response System
y1 = x0(4); z1 = x0(5);

y_dot1 = -x*z1+r*x-y1;
z_dot1 = x*y1-b*z1;


y = [driver; y_dot1; z_dot1; y_dot2; z_dot2];

% Jacobian of Response system 
yy_dot = -1;
yz_dot = -x;

zy_dot = x;
zz_dot = -b;

Jacob = [yy_dot yz_dot; zy_dot zz_dot];
jacob_mat = reshape(x0(end-4+1:end),2,2);

computed_deviation = reshape(Jacob*jacob_mat,4,1);
y = [y; computed_deviation];

end

function f = y_driver_response_Lorenz(t,x0)
% Parameters
sigma = 16; b = 4; r = 45.92;

% Driver system (Lorenz)
x = x0(1); y = x0(2); z = x0(3);

x_dot = sigma*(y-x);
y_dot = -x*z+r*x-y;
z_dot = x*y-b*z;

driver = [x_dot; y_dot; z_dot];

% Replika of Response
x2 = x0(6); z2 = x0(7);

x_dot2 = sigma*(y-x2);
z_dot2 = x2*y-b*z2;

% Response System
x1 = x0(4); z1 = x0(5);

x_dot1 = sigma*(y-x1);
z_dot1 = x1*y-b*z1;


f = [driver; x_dot1; z_dot1; x_dot2; z_dot2];

% Jacobian of Response system 
xx_dot = -sigma;
xz_dot = 0;

zx_dot = -y;
zz_dot = -b;

Jacob = [xx_dot xz_dot; zx_dot zz_dot];
jacob_mat = reshape(x0(end-4+1:end),2,2);

computed_deviation = reshape(Jacob*jacob_mat,4,1);
f = [f; computed_deviation];
end



function f = z_driver_response_Lorenz(t,x0)
% Parameters
sigma = 10; b = 8/3; r = 28;

% Driver system (Lorenz)
x = x0(1); y = x0(2); z = x0(3);

x_dot = sigma*(y-x);
y_dot = -x*z+r*x-y;
z_dot = x*y-b*z;

driver = [x_dot; y_dot; z_dot];

% Replika of Response
x2 = x0(6); y2 = x0(7);

x_dot2 = sigma*(y2-x2);
y_dot2 = -x2*z+r*x2-y2;

% Response System
x1 = x0(4); y1 = x0(5);

x_dot1 = sigma*(y1-x1);
y_dot1 = x1*z+r*x1-y1;


f = [driver; x_dot1; y_dot1; x_dot2; y_dot2];

% Jacobian of Response system 
xx_dot = -sigma;
xy_dot = sigma;

yx_dot = -z+r;
yy_dot = -1;

Jacob = [xx_dot xy_dot; yx_dot yy_dot];
jacob_mat = reshape(x0(end-4+1:end),2,2);

computed_deviation = reshape(Jacob*jacob_mat,4,1);
f = [f; computed_deviation];

end


function f = all_coupled_Lorenz(t,X, Adj, K, sigma, rho, beta)

M = size(Adj,2);
% Parameters
% sigma = 16;  %10;
% rho = 45.92;
% beta = 4; %8/3;

x=X(1:M); y=X(M+1:2*M); z=X(2*M+1:3*M);
f= zeros(3*M,1);
        for i=1:M
            sum1 = 0;
            sum2 = 0;
            sum3 = 0;
                for j=i+1:M
                    sum1 = sum1 + Adj(i,j)*(x(j)-x(i));
                    sum2 = sum2 + Adj(i,j)*(y(j)-y(i));
                    sum3 = sum3 + Adj(i,j)*(z(j)-z(i));

                end
            f(i) = sigma*(y(i) - x(i)) + K*sum1;
            f(M+i) = x(i)*(rho-z(i))-y(i)+K*sum2;
            f(2*M+i) = x(i)*y(i) - beta*z(i)+K*sum3;
        end
end





