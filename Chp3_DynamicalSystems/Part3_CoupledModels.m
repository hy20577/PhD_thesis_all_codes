
%% Part3: Coupled Systems

clc; clear; close;

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp1 = strjoin(splittedstr(1:end-1), '/');
addpath(genpath(sprintf('%s/Functions', pp1)))


fig_folder = sprintf('%s/Figures',current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end
%% 9. Network Figures
Adj1 = readmatrix(sprintf('%s/Functions/12nodes.txt', pp1));
Adj2 = readmatrix(sprintf('%s/Functions/16nodes.txt', pp1));

fig = figure;
sb1 = subplot(1,2,1);
plot(graph(Adj1), 'MarkerSize',12, 'NodeColor','r', 'EdgeColor','k', 'LineWidth', 2, 'Layout','circle', 'NodeFontSize',20)
ax = gca;
ax.Box = 'off';
% xlim([-1.2,1.2]); ylim([-1.2,1.2])

ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
axis equal
title('(a)', 'FontSize',30, 'Interpreter','latex')
pbaspect([1 1 1])

sb2 = subplot(1,2,2);
plot(graph(Adj2), 'MarkerSize',12, 'NodeColor','r', 'EdgeColor','k', 'LineWidth', 2, 'Layout', 'circle', 'NodeFontSize',20)
ax = gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
pbaspect([1 1 1])
title('(b)', 'FontSize',30, 'Interpreter','latex')

%%

saveas(gcf, sprintf('%s/NetworkTopology.fig', fig_folder))
saveas(gcf, sprintf('%s/NetworkTopology', fig_folder), 'epsc')

%% 10. Bifurcation and MLE for Logistic & Circle

k=1e4;
a=linspace(1,4,k);
t=1e3;
B = zeros(100,k);
Adj = 1;
niter= 1e4;
LYA = zeros(k,1);

parfor j=1:k
    x0=rand(1,1);
    [Orbits, Lya_overtime] = Log_generator(Adj, 0, x0, a(j), niter, 2, 1, 0, 0); 
    B(:,j) = Orbits(end-99:end);
    LYA(j) = Lya_overtime(end);
end
%
figure
subplot(2,1,1)
for m=1:k
plot(repmat(a(m),t-900,1),5*B(:,m),'k.')
hold on
plot(a(m), LYA(m), 'r.')
end
ylim([-10,5])
yline(0, 'LineWidth',2)

%
ax = gca;
ax.YTickLabel(end) = {'1'};
ax.FontSize = 30;
xlabel('$r$', 'interpreter', 'latex', 'FontSize',30)
title('(a)', 'interpreter', 'latex', 'FontSize',30)
%
% saveas(gcf, sprintf('%s/LogisticBifurcationMLE.fig', fig_folder))
% saveas(gcf, sprintf('%s/LogisticBifurcationMLE', fig_folder), 'epsc')
%  Bifurcation and MLE for Circle
%
k=3e3;
ac=linspace(0,1,k);
t=1e3;
Bc = zeros(100,k);
Adj = 1;
niter= 1e4;
LYA2 = zeros(k,1);

parfor j=1:k
    x0=rand(1,1);
    [Orbits, Lya_overtime] = Circle_generator(Adj, 0, x0, ac(j), niter, 2, 1, 0, 0); 
    Bc(:,j) = Orbits(end-99:end);
    LYA2(j) = Lya_overtime(end);
end
%
subplot(2,1,2)
for m=1:k
plot(repmat(ac(m),t-900,1),Bc(:,m),'k.')
hold on
plot(ac(m), LYA2(m), 'r.')
end
yline(0, 'Linewidth',2)
ax2 = gca;
ax2.FontSize = 30;
%%
ax2.Title.String = '(b)';
ax2.Title.FontSize = 30;
ax2.Title.Interpreter ='latex';
ax2.XLabel.String = '$r$';
ax2.XLabel.FontSize = 30;
ax2.XLabel.Interpreter ='latex';

%%
saveas(gcf, sprintf('%s/LogCircBifurcationMLE.fig', fig_folder))
saveas(gcf, sprintf('%s/LogCircBifurcationMLE', fig_folder), 'epsc')


%%  %% %%%%%%%%%%%% 11. Coupled HR  %%%%%%%%%%%

% Uncoupled

t0 = 0;
tf = 1e3;
k = 1e4;
dt = (tf-t0)/k;
A = readmatrix(sprintf('%s/Functions/12nodes.txt', pp1));

%Initial Conditions
M = length(A(1,:));
etap = 0.5*rand(1,M);
etaq = 0.5*rand(1,M);
etan = 0.5*rand(1,M);

p_init = zeros(1,M);
q_init = zeros(1,M);
n_init = zeros(1,M);

for i=1:M
        p_init(i)   = -1.30784489 + etap(i);
        q_init(i)   = -7.32183132 + etaq(i);      %% Setting initial conditions
        n_init(i)   =  3.35299859 + etan(i);
end

x0 = [p_init q_init n_init];
coupling_strength = 0;
time_window = 5e2;
std_threshold = 1e-3;
max_discard= 1e4;

[tvalues, LYA, Orbits] = HR_generator(x0, t0, tf, k, coupling_strength, A, time_window, std_threshold, max_discard);

%
LYA_uncoupled = LYA(end,:);

inf = readmatrix(sprintf('HR_Info/info_%.2f.txt', coupling_strength));
disc = inf(1,2);
%

%
% 0.1 Coupled

coupling_strength = 0.1;

[tvalues_01, LYA_01, Orbits_01] = HR_generator(x0, t0, tf, k, coupling_strength, A, time_window, std_threshold, max_discard);

%
LYA_01_coupled = LYA_01(end,:);
inf = readmatrix(sprintf('HR_Info/info_%.2f.txt', coupling_strength));
disc1 = inf(1,2);
%%

figure
subplot(3,2,1)
plot3(Orbits(:,1),Orbits(:,M+1),Orbits(:,2*M+1))
xlabel('p', 'interpreter', 'latex', 'Fontsize',20); 
ylabel('q', 'interpreter', 'latex', 'Fontsize',20); 
zlabel('n', 'interpreter', 'latex', 'Fontsize',20)
title('(a)', 'FontSize',20, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30;
xticks([-2,2]); yticks([-10, 0]); zticks([2.5,3.5])

xlim([-2,2]); ylim([-12,1]); zlim([2.5,3.8])

subplot(3,2,3)
plot(tvalues(disc+1:end), Orbits(:,1), 'LineWidth',2)
hold on
plot(tvalues(disc+1:end), Orbits(:,M+1), 'LineWidth',2)
plot(tvalues(disc+1:end), Orbits(:,2*M+1), 'LineWidth',2)
% legend('p', 'q', 'n', 'interpreter', 'latex', 'Fontsize', 12, 'Location', 'bestoutside')
% ylabel('Magnitude')
xlim([0 16e2]); ylim([-12, 5]);
title('(b)', 'FontSize',20, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30;
ax.XTick = [];
ax.YTick = [ax.YTick(1), ax.YTick(end)];
ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};

subplot(3,2,5)
semilogy(0.1*(1:size(LYA,1)), LYA(:,1), 'LineWidth',2)
xline(0.1*disc, 'r--', 'LineWidth',2)
xlabel('Time', 'FontSize',20, 'interpreter', 'latex'); 
ylabel('MLE', 'FontSize',20, 'interpreter', 'latex'); xlim([0 16e2]);
title('(c)', 'FontSize',20, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30;


subplot(3,2,2)
plot3(Orbits_01(:,1),Orbits_01(:,M+1),Orbits_01(:,2*M+1))
xlabel('p', 'interpreter', 'latex', 'Fontsize',20); 
ylabel('q', 'interpreter', 'latex', 'Fontsize',20); 
zlabel('n', 'interpreter', 'latex', 'Fontsize',20)
title('(d)', 'FontSize',20, 'Interpreter','latex');
xlim([-2,2]); ylim([-12,1]); zlim([2.5,3.8])
ax = gca;
ax.FontSize = 30;
ax.XTick = [ax.XTick(1), ax.XTick(end)];
ax.XTickLabel = {sprintf('%s', ax.XTickLabel{1}),sprintf('%s', ax.XTickLabel{end}) };
ax.YTick = [ax.YTick(1), ax.YTick(end)];
ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};
ax.ZTick = [ax.ZTick(1), ax.ZTick(end)];
ax.ZTickLabel = {sprintf('%s', ax.ZTickLabel{1}),sprintf('%s', ax.ZTickLabel{end})};


subplot(3,2,4)
plot(tvalues_01(disc1+1:end), Orbits_01(:,1), 'LineWidth',2)
hold on
plot(tvalues_01(disc1+1:end), Orbits_01(:,M+1), 'LineWidth',2)
plot(tvalues_01(disc1+1:end), Orbits_01(:,2*M+1), 'LineWidth',2)
legend('p', 'q', 'n', 'interpreter', 'latex', 'Fontsize', 20, 'Location', 'bestoutside')
% ylabel('Magnitude')
xlim([0 14e2]); ylim([-12, 5]);
title('(e)', 'FontSize',20, 'Interpreter','latex')
ax = gca;
ax.XTick = [];
ax.FontSize = 30;
ax.YTick = [ax.YTick(1), ax.YTick(end)];
ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};

subplot(3,2,6)
semilogy(0.1*(1:size(LYA_01,1)), LYA_01(:,1), 'LineWidth',2)
xline(0.1*disc1, 'r--', 'LineWidth',2); xlim([0 14e2]);
xlabel('Time', 'FontSize',20, 'interpreter', 'latex')
title('(f)', 'FontSize',20, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/HR_coupled_uncoupled.fig', fig_folder))
saveas(gcf, sprintf('%s/HR_coupled_uncoupled', fig_folder), 'epsc')

%% 12. Coupled Lorenz
Adj = [0 1;
       1 0];
M = size(Adj,2);
t0=0; tf=2e5; k=4e5; 
sigma=10; beta = 8/3; rho = 28;
coupling_strength = linspace(0,1.5,16);
x_init = 0.1 + 0.5*rand(1,M);
y_init = -0.19 + 0.5*rand(1,M);
z_init = -0.27 + 0.5*rand(1,M);
xyz0 = [x_init, y_init, z_init]; 
time_window=5e2;
std_threshold=1e-2;
max_discard=1e4;

parfor i=1:length(coupling_strength)
    [tvalues, LYA_Lorenz, Orbits_Lorenz] = Lorenz_generator(xyz0, t0, tf, k, coupling_strength(i), sigma, rho, beta, Adj, time_window, std_threshold, max_discard);
    MLE(i) = max(LYA_Lorenz(end,:));
end
%%
coupling_strength = linspace(0,1.5,16);
MLE = zeros(16,1);
for i=1:16
m = readmatrix(sprintf('Lorenz_Info/MLE_%.2f.txt',round(coupling_strength(i),2)));
MLE(i) = m(end);
end
Lorenz_coup_MLE = [coupling_strength', MLE];
writematrix(Lorenz_coup_MLE, 'Lorenz_coup_MLE.txt')

%%
t0=0; tf=2e3; k=2e5; 
[tvalues, LYA_Lorenz, Orbits_Lorenz] = Lorenz_generator(xyz0, t0, tf, k, 0, sigma, rho, beta, Adj, time_window, std_threshold, max_discard);


%%

figure
subplot(2,2,4)
plot(coupling_strength, MLE, '-', 'Linewidth', 2)
hold on
plot(coupling_strength, MLE, 'r*', 'Linewidth', 2)
ylim([0.5,1])
% idx = find(MLE == max(MLE));
% hold on
% plot(coupling_strength(idx), MLE(idx), 'o')
xlabel('K',  'interpreter', 'latex', 'Fontsize',20)
ylabel('MLE',  'interpreter', 'latex', 'Fontsize',20)
title('(c)', 'interpreter', 'latex', 'FontSize',20)
ax= gca;
ax.FontSize = 30;
ax.XTick = [ax.XTick(1), ax.XTick(end)];
ax.XTickLabel = {sprintf('%s', ax.XTickLabel{1}),sprintf('%s', ax.XTickLabel{end}) };
ax.YTick = [ax.YTick(1), ax.YTick(end)];
ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};

subplot(2,2,2)
loglog((1:size(LYA_Lorenz,1))*0.1, LYA_Lorenz(:,1), 'LineWidth',2)
xline(0.1*(size(LYA_Lorenz,1)-k), 'r--', 'LineWidth',2)
xlabel('Time', 'Interpreter','latex', 'FontSize',20);
ylabel('MLE', 'Interpreter','latex', 'FontSize',20)
title('(b)', 'interpreter', 'latex', 'FontSize',20)
ax= gca;
ax.FontSize = 30;
ax.XTick = [ax.XTick(1), ax.XTick(end)];
ax.XTickLabel = {sprintf('%s', ax.XTickLabel{1}),sprintf('%s', ax.XTickLabel{end}) };
% ax.YTick = [ax.YTick(1), ax.YTick(end)];
% ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};

subplot(2,2,[1,3])
plot3(Orbits_Lorenz(end-1e5:end,1),Orbits_Lorenz(end-1e5:end,M+1),Orbits_Lorenz(end-1e5:end,2*M+1))
xlabel('x', 'interpreter', 'latex', 'Fontsize',20); ylabel('y', 'interpreter', 'latex', 'Fontsize',20); zlabel('z', 'interpreter', 'latex', 'Fontsize',20);
title('(a)', 'interpreter', 'latex', 'FontSize',20)
ax= gca;
ax.FontSize = 30;

ax.XTick = [ax.XTick(1), ax.XTick(end)];
ax.XTickLabel = {sprintf('%s', ax.XTickLabel{1}),sprintf('%s', ax.XTickLabel{end}) };
ax.YTick = [ax.YTick(1), ax.YTick(end)];
ax.YTickLabel = {sprintf('%s', ax.YTickLabel{1}),sprintf('%s', ax.YTickLabel{end})};
ax.ZTick = [ax.ZTick(1), ax.ZTick(end)];
ax.ZTickLabel = {sprintf('%s', ax.ZTickLabel{1}),sprintf('%s', ax.ZTickLabel{end})};


%%
saveas(gcf, sprintf('%s/LorenzMLE_trajectory.fig',fig_folder));
saveas(gcf, sprintf('%s/LorenzMLE_trajectory',fig_folder), 'epsc');

%% %% 13. Deterministic Kuramoto Generator 

% From the written .txt file
n= 32;
coupling_strength = linspace(2,10,n);
folder = 'Kuramoto_Info';
MLEs = zeros(n,1);

for i=1:n
    info = readmatrix(sprintf('%s/info_%.2f.txt',folder, coupling_strength(i)));
    MLEs(i) = info(2,2); 
end

%%

% figure
% subplot(2,3,[1 4])
% plot(coupling_strength, MLEs, 'o--')  % For 2 connected nodes.
% ylabel('MLE', 'interpreter', 'latex', 'Fontsize',12)
% xlabel('K', 'interpreter', 'latex', 'Fontsize',12)

%

Adj = readmatrix(sprintf('%s/Functions/12nodes.txt', pp1));
M = size(Adj,2);
n= 32;
coupling_strength = linspace(2,10,n);

% % % % % % % %  Data generation over coupling strength % % % % % % % % %

x0 = rand(M,1); % randomly chosen M initial conditions
t0 = 0;
tf = 5e3;
k = 5e4;
w = rand(1,M)*2*pi-pi;   %internal frequencies has 0 mean, [-pi, pi].
time_window=5e2;
std_threshold=1e-2;
max_discard=1e4;

MLEs = zeros(length(coupling_strength),1);

parfor i=1:n
    [~, LYA, ~] = Kuramoto_generator(x0, t0, tf, k, coupling_strength(i), Adj, w, time_window, std_threshold, max_discard);
    MLEs(i) = max(LYA(end,:));
end
%%
id = find(MLEs == max(MLEs));
MLEs(id)
coupling_strength(id)

coup = coupling_strength(id);

% % % % % % % %  Data generation for the maximum LYA  % % % % % % % % %

x0 = rand(M,1); % randomly chosen M initial conditions
t0 = 0;
tf = 5e3;
k = 5e4;
w = rand(1,M)*2*pi-pi;   %internal frequencies has 0 mean, [-pi, pi].
time_window=5e2;
std_threshold=1e-2;
max_discard=1e4;

[tvalues, LYA, Orbits] = Kuramoto_generator(x0, t0, tf, k, coup, Adj, w, time_window, std_threshold, max_discard);

inst_freq = zeros(size(Orbits,1)-1, size(Orbits,2));
% 
for i=1:k-1
    inst_freq(i,:) = (Orbits(i+1,:)- Orbits(i, :))/(tvalues(2)- tvalues(1));
end

info = readmatrix(sprintf('Kuramoto_Info/info_%.2f.txt', coup));
disc = info(1,2); 

figure
subplot(2,2,1)
semilogy(coupling_strength, MLEs, '*-', 'LineWidth',2)  % For 2 connected nodes.
ylabel('MLE', 'interpreter', 'latex', 'Fontsize',20)
xlabel('K', 'interpreter', 'latex', 'Fontsize',20)
hold on
plot(coupling_strength(id), MLEs(id), 'ro', 'LineWidth',8)
title('(a)', 'interpreter', 'latex', 'Fontsize',20)
ax= gca;
ax.FontSize = 30;

subplot(2,2,2)
plot(tvalues(disc+1:end-48e3), Orbits(1:end-48e3,1), 'LineWidth',2)
hold on
plot(tvalues(disc+1:end-48e3), Orbits(1:end-48e3,2), 'LineWidth',2, 'Color','#eb8c34')
plot(tvalues(disc+1:end-48e3), Orbits(1:end-48e3,3), 'LineWidth',2)
plot(tvalues(disc+1:end-48e3), Orbits(1:end-48e3,4), 'LineWidth',2)
plot(tvalues(disc+1:end-48e3), Orbits(1:end-48e3,5), 'LineWidth',2, 'Color', '#3480eb')
xlim([0 3e2])
xlabel('Time', 'interpreter', 'latex', 'Fontsize',20)
title('(c)', 'interpreter', 'latex', 'Fontsize',20)
legend('1', '2', '3', '4', '5', 'Location', 'bestoutside')
ax= gca;
ax.FontSize = 30;

subplot(2,2,4)
plot(tvalues(disc+2:end-48e3), inst_freq(1:end-48e3,2), 'Color','#eb8c34')
hold on
plot(tvalues(disc+2:end-48e3), inst_freq(1:end-48e3,5), 'Color', '#3480eb')
xlim([0 3e2])
xlabel('Time', 'interpreter', 'latex', 'Fontsize',20)
title('(d)', 'interpreter', 'latex', 'Fontsize',20)
ax= gca;
ax.FontSize = 30;

subplot(2,2,3)
semilogy(tvalues, LYA(:,1))
xline(disc*0.1, 'r--')
title('(b)', 'interpreter', 'latex', 'Fontsize',20)
xlabel('Time', 'interpreter', 'latex', 'Fontsize',20)
ylabel('MLE', 'interpreter', 'latex', 'Fontsize',20)
ax= gca;
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/DeterministicKuramoto.fig', fig_folder))
saveas(gcf, sprintf('%s/DeterministicKuramoto', fig_folder), 'epsc')

%%
%% 14. Stochastic Kuramoto Data Generalization

% global A D w

% A= Adj;
nVar = length(A(1,:));

w = 2*pi*rand(1,nVar)-pi;  %% Internal frequency of the oscilators.
dt=0.1;
t0=0;
tf=5e4;
niter = (tf-t0)/dt;
X=  2*pi*rand(nVar,1);

D = 10;
 
obj = sde(@drift, @diff, 'StartState', X);
%%
[iteration, tvals] = simulate(obj, niter, 'DeltaTime', dt);


Instfreq = zeros(length(tvals)-1,nVar);
for i=1:length(tvals)-1
    Instfreq(i,:) = (iteration(i+1,:) - iteration(i,:))/dt;
end


figure
subplot(2,1,1)
plot(tvals(1:5e3), iteration(1:5e3,1), 'LineWidth',2, 'Color', '#a8232a')
hold on
plot(tvals(1:5e3), iteration(1:5e3,2), 'LineWidth',2, 'Color','black')
plot(tvals(1:5e3), iteration(1:5e3,3), 'LineWidth',2)
plot(tvals(1:5e3), iteration(1:5e3,4), 'LineWidth',2)
plot(tvals(1:5e3), iteration(1:5e3,5), 'LineWidth',2)
xlabel('Time', 'interpreter', 'latex', 'FontSize',20)
title('(a)', 'interpreter', 'latex', 'FontSize',20)
legend('1', '2', '3', '4', '5', 'Location', 'bestoutside')
ax= gca;
ax.FontSize = 30; 

subplot(2,1,2)
plot(tvals(2:5e3), Instfreq(1:5e3-1,2), 'Color','black')
hold on
plot(tvals(2:5e3), Instfreq(1:5e3-1,1), 'Color', '#a8232a')

xlabel('Time', 'interpreter', 'latex', 'FontSize',20)
title('(b)', 'interpreter', 'latex', 'FontSize',20)
ax= gca;
ax.FontSize =30; 


%%
saveas(gcf, sprintf('%s/StochasticKuramoto.fig', fig_folder))
saveas(gcf, sprintf('%s/StochasticKuramoto', fig_folder), 'epsc')




%%%%%%%%%%%  FUNCTIONS  %%%%%%

function dr = drift(t,x, A, w)
% A is the adjaceny matrix
% theta is phases of oscilators
% K is the coupling strength


% global A w

K=10; 
N = size(A,2);
dr = zeros(N,1);

    for i=1:N
        sum1 =0; 
        for j=1:N
            sum1 = sum1+ A(i,j) * sin(x(j)-x(i));  
        end
        dr(i) = w(i) + (K/N)*sum1; 
    end
end
% 
function diffusion = diff(t,x, A, D)
% global A D
N = length(A(1,:));
diffusion= D*eye(N);
end





