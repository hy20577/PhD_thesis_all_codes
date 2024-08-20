
close all; clear; clc;
%
current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp = strjoin(splittedstr(1:end-1), '/'); 

addpath(genpath(sprintf(pp, '/Functions')))

fig_folder = sprintf('%s/Figures', current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end

%% Ring Lattice

A = ringAdjMat(16,4);
nedge = sum(sum(A));
N = size(A,1);
dist_mat = distances(graph(A));
dist_array = zeros(N,1);

counter = 0;
for i=1:N
    for j= i+1:N
        counter = counter +1;
        dist_array(counter) = dist_mat(i,j);
    end
end

ASPL = mean(dist_array);
diameter = max(max(dist_array));

figure
subplot(2,2,1)
h = plot(graph(A), 'MarkerSize',20, 'NodeColor','#f5a742', 'NodeFontSize',15, 'LineWidth',2, 'EdgeColor','#7a7c80', 'Layout','circle');
pbaspect([1 1 1])
 
text(h.XData, h.YData,  h.NodeLabel, 'VerticalAlignment','middle',...
 'HorizontalAlignment', 'center', 'FontSize', 30)
h.NodeLabel = {}; 
ax=gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis

title('(a)', 'FontSize',30, 'Interpreter','latex')

subplot(2,2,3)
histogram(sum(A), 'Normalization','probability', 'FaceColor', '#f5a742', 'EdgeColor','k', 'LineWidth',2)
ylim([0,1])
xlim([0,10])
pbaspect([1 1 1])
title('(b)', 'FontSize',30, 'Interpreter','latex')
xlabel('$k$', 'FontSize',30, 'Interpreter','latex')
ylabel('$p(k)$', 'FontSize',30, 'Interpreter','latex')
ax=gca;
ax.Box = 'off';
ax.FontSize = 30;
yticks([0, 0.5, 1])

N = 10:1:40;
average_degree = zeros(length(N),1);
shortest = zeros(length(N),1); 
counter = 0;
for i=1:length(N)
A = ringAdjMat(N(i), 3);
counter = counter +1; 
average_degree(counter) = sum(sum(distances(graph(A))))/(N(i)*(N(i)-1));
shortest(counter) = max(max(distances(graph(A)))); 
end

subplot(2,2,[2,4])
plot(N, average_degree, 'k', 'Linewidth', 2)
hold on
plot(N, shortest, 'r', 'LineWidth',2)
legend('$\overline{d}$', '$R$', 'interpreter', 'latex')
title('(c)', 'FontSize',30, 'Interpreter','latex')
xlabel('$N$', 'FontSize',30, 'Interpreter','latex')
ax = gca; 
ax.FontSize = 30; 
%%
saveas(gcf, 'regularlattice.fig')
saveas(gcf, 'regularlattice', 'epsc')

%% NETWORK GRAPHS AND DEGREE DISTRIBUTIONS 

ER = readmatrix('Adjacency/Erdos Renyi.txt');
SW = readmatrix('Adjacency/Small World.txt');
SF = readmatrix('Adjacency/Scale Free.txt'); 
%%
% number of edges

v_ER = sum(sum(ER))/2;
v_SW = sum(sum(SW))/2;
v_SF = sum(sum(SF))/2;

% density
N = 40;
allpairs = 40*(39)/2;
d_ER = v_ER/allpairs;
d_SW = v_SW/allpairs;
d_SF = v_SF/allpairs;

dist_matrix_ER = distances(graph(ER));
dist_matrix_SW = distances(graph(SW));
dist_matrix_SF = distances(graph(SF)); 

diameter_ER = max(max(dist_matrix_ER));
diameter_SW = max(max(dist_matrix_SW));
diameter_SF = max(max(dist_matrix_SF));

dist_array_ER = zeros(allpairs,1);
dist_array_SW = zeros(allpairs,1);
dist_array_SF = zeros(allpairs,1);

counter = 0;
for i=1:40
    for j=i+1:40
        counter = counter +1;
        dist_array_ER(counter) = dist_matrix_ER(i,j);
        dist_array_SW(counter) = dist_matrix_SW(i,j);
        dist_array_SF(counter) = dist_matrix_SF(i,j);
    end
end

ASPL_ER = mean(dist_array_ER);
ASPL_SW = mean(dist_array_SW);
ASPL_SF = mean(dist_array_SF);


%%
figure
subplot(1,2,1)
plot(graph(ER), 'MarkerSize',20, 'NodeColor','#f5a742', 'NodeFontSize',20, 'LineWidth',2, 'EdgeColor','#7a7c80')
pbaspect([1 1 1])
ax=gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
% axis equal 
title('(a)', 'FontSize',30, 'Interpreter','latex')

subplot(1,2,2)
histogram(sum(ER), 'FaceColor','#f5a742', 'LineWidth',2)
pbaspect([1 1 1])
ax= gca; 
ax.FontSize = 30; 
title('(b)', 'FontSize',30, 'Interpreter','latex')
xlabel('Node Degree', 'FontSize',30, 'Interpreter','latex')
%%
saveas(gca, sprintf('%s/ER_network.fig', fig_folder))
saveas(gca, sprintf('%s/ER_network', fig_folder), 'epsc')

%%
figure
subplot(1,2,1)
plot(graph(SW), 'MarkerSize',20, 'NodeColor','#f5a742', 'NodeFontSize',20, 'LineWidth',2, 'EdgeColor','#7a7c80')
pbaspect([1 1 1])
ax=gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
title('(a)', 'FontSize',30, 'Interpreter','latex')

subplot(1,2,2)
histogram(sum(SW),'Normalization', 'probability', 'FaceColor','#f5a742', 'LineWidth',2)
pbaspect([1 1 1])
ax= gca; 
ax.Box = 'off';
ax.FontSize = 30; 
title('(b)', 'FontSize',30, 'Interpreter','latex')
xlabel('Node Degree', 'FontSize',30, 'Interpreter','latex')
ylim([0,1])
xlabel('$k$', 'FontSize',30, 'Interpreter','latex')
ylabel('$p(k)$', 'FontSize',30, 'Interpreter','latex')
%%
saveas(gca, sprintf('%s/SW_network.fig', fig_folder))
saveas(gca, sprintf('%s/SW_network', fig_folder), 'epsc')

%%
figure
subplot(1,2,1)
plot(graph(SF), 'MarkerSize',20, 'NodeColor','#f5a742', 'NodeFontSize',20, 'LineWidth',2, 'EdgeColor','#7a7c80')
pbaspect([1 1 1])
ax=gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off'; % Hide the x-axis
ax.YAxis.Visible = 'off'; % Hide the y-axis
title('(a)', 'FontSize',30, 'Interpreter','latex')
ax.Box = 'off';

subplot(1,2,2)
histogram(sum(SF),'Normalization', 'probability', 'FaceColor','#f5a742', 'LineWidth',2)
pbaspect([1 1 1])
ax= gca; 
ax.Box = 'off';
ax.FontSize = 30; 
title('(b)', 'FontSize',30, 'Interpreter','latex')
%xlabel('Node Degree', 'FontSize',30, 'Interpreter','latex')
ylim([0,1])
xlabel('$k$', 'FontSize',30, 'Interpreter','latex')
ylabel('$p(k)$', 'FontSize',30, 'Interpreter','latex')
hold on
plot(1:0.1:20, (1:0.1:20).^-1, 'LineWidth',2)
%%
saveas(gca, sprintf('%s/SF_network.fig', fig_folder))
saveas(gca, sprintf('%s/SF_network', fig_folder), 'epsc')

















