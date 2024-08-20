
close, clear, clc

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp = strjoin(splittedstr(1:end-1), '/'); 
addpath(genpath(sprintf(pp, '/Functions')))
% 

fig_folder = sprintf('%s/Figures', current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end

%% %% %% %% %% Phase sync. small-coupled rossler,

%% Rossler phase synchronization

Adj  = [0 1 0 0 0;
        1 0 1 1 0;
        0 1 0 0 0;
        0 1 0 0 1;
        0 0 0 1 0];
M = size(Adj,2);
x0 = rand(3*M,1);

coupling_strength = 0.02;
a = 0.15; f= 0.2; c= 10;
dw = 0.015; % frequency mismatch
w = linspace(1-dw, 1+dw, 5); % 1+2*rand(M,1)*dw-dw;

t0 = 0;
tf = 1e4;
k=1e5;
dt = (tf-t0)/(k-1);

time_window = 5e2;
std_threshold = 1e-3;
max_discard = 1e4;


[tvalues, LYA, Orbits] = rossler_generator( x0, t0, tf, k, coupling_strength, a, f, c, w, Adj, time_window, std_threshold, max_discard);
inst_phases = hilbert_phases(Orbits);

%%
writematrix(Orbits, 'Rossler_Info/Orbits_rossler.txt')
writematrix(tvalues, 'Rossler_Info/tvals_rossler.txt')
writematrix(LYA, 'Rossler_Info/LYA.txt')
writematrix(inst_phases, 'Rossler_Info/inst_phases.txt')

%%
% Instantaneous phases
coupling_strength = 0.02;
Orbits = readmatrix('Rossler_Info/Orbits_rossler.txt'); 
inst_phases = readmatrix('Rossler_Info/inst_phases.txt');
tvalues = readmatrix('Rossler_Info/tvals_rossler.txt');
info =readtable(sprintf('Rossler_Info/info_%.3f.txt', coupling_strength));
discarded = info{1,2};


%

figure
subplot(1,2,1)
plot(tvalues(discarded+1:end), abs(inst_phases(:,1)-inst_phases(:,2)), 'LineWidth',2)
hold on
plot(tvalues(discarded+1:end), abs(inst_phases(:,1)-inst_phases(:,3)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,2)-inst_phases(:,3)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,4)-inst_phases(:,5)), 'LineWidth',2)
ylim([0,150])
xlabel('Time', 'interpreter', 'latex', 'Fontsize', 30)
ylabel('$\psi_{(1,1)}$', 'interpreter', 'latex', 'Fontsize', 30)
legend('1-2', '1-3', '2-4', '4-5','interpreter', 'latex', 'Fontsize', 30, 'location', 'best')
xticks([0, 10000]);
xticklabels({'0', '10^{4}'})
title('(a)', 'interpreter', 'latex', 'FontSize',30)
ax= gca;
ax.FontSize = 30;

subplot(1,2,2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,1)-inst_phases(:,4)), 'LineWidth',2)
hold on
plot(tvalues(discarded+1:end), abs(inst_phases(:,1)-inst_phases(:,5)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,2)-inst_phases(:,4)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,2)-inst_phases(:,5)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,3)-inst_phases(:,4)), 'LineWidth',2)
plot(tvalues(discarded+1:end), abs(inst_phases(:,3)-inst_phases(:,5)), 'LineWidth',2)
ylim([0,150])
legend('1-4', '1-5', '2-3', '2-5', '3-4', '3-5','interpreter', 'latex', 'Fontsize', 30, 'location', 'best')
title('(b)', 'interpreter', 'latex', 'FontSize',30)
xlabel('Time', 'interpreter', 'latex', 'Fontsize', 30)
ylabel('$\psi_{(1,1)}$', 'interpreter', 'latex', 'Fontsize', 30)
xticks([0, 10000]);
xticklabels({'0', '10^{4}'})
ax= gca;
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/phase_sync_rossler.fig', fig_folder))
saveas(gcf, sprintf('%s/phase_sync_rossler', fig_folder), 'epsc')
%%  Magnitude Order and Pairwise Order  (Logistic Map of 5-nodes network)

Lya_overtime = readmatrix('Logistic_Info/MLE_0.20.txt');
Orbits = readmatrix('Logistic_Info/Orbits.txt'); 
info = readmatrix('Logistic_Info/info_0.20.txt');
discarded = info(1,2);
M = size(Orbits,2); 
% global_order_logistic = readmatrix('Logistic_Info/global_order_logistic.txt'); 
% pcs = readmatrix('Logistic_Info/pws_order_logistic.txt'); 
%
%%
pcs = zeros(M,M);
for i=1:M
    for j=i+1:M
        R = globalorder_map(Orbits(:,[i,j]));
        pcs(i,j) = R(end);
    end
end
%
pcs = pcs + pcs';
%%
all_globalorder = globalorder_map(Orbits);

%%
writematrix(all_globalorder, 'Logistic_Info/global_order_logistic.txt')
writematrix(pcs, 'Logistic_Info/pws_order_logistic.txt')
%
% pcs2 = NaN(M,M);
% for i=1:M
%     for j=i+1:M
%         R2 = globalorder_map(Orbits(1:5e4,[i,j]));
%         pcs2(i,j) = R2(end);
%     end
% end
% 
% %%
% 
% [c1, p1] = corr(Orbits);

%%
for i=1:5; pcs(i,i) = NaN; end

figure
subplot(1,2, 1)
plot(all_globalorder, 'LineWidth',2) % global order
title('(a)', 'interpreter', 'latex', 'Fontsize', 30);
xlabel('Length', 'interpreter', 'latex', 'Fontsize', 30)
ylabel('R', 'interpreter', 'latex', 'Fontsize', 30)
xticks([0, 5e4, 1e5]);
yticks([0, 0.2, 0.4])
ax = gca;
ax.FontSize = 30;
pbaspect([1 1 1])

subplot(1,2,2)
h1 = imagesc(pcs);
clbar = colorbar();
% title('(b)', 'interpreter', 'latex', 'FontSize',30)
pbaspect([1 1 1])
% ax.Title = '(b)';
% ax.Title.Interpreter = 'latex'; 
clbar.Ticks = [0.5, 0.6, 0.7, 0.8];
ax = gca;
ax.FontSize = 30;
title('(b)', 'interpreter', 'latex', 'FontSize',30)
% subplot(2,2,3)
% h1 = imagesc(c1);
% colorbar();

% pbaspect([1 1 1])
% % ax.Title = '(d)';
% % ax.Title.Interpreter = 'latex'; 
% ax = gca;
% ax.FontSize = 30;
% 
% subplot(2,2,4)
% h1 = imagesc(p1);
% colorbar();
% 
% % title('(b)', 'interpreter', 'latex', 'FontSize',30)
% pbaspect([1 1 1])
% % ax.Title = '(d)';
% % ax.Title.Interpreter = 'latex'; 
% ax = gca;
% ax.FontSize = 30;
%%
% set(gcf,'color','w'); % white background
% set(gcf, 'Position', get(0,'ScreenSize'))
saveas(gcf, sprintf('%s/cs_magnitudeOrder.fig', fig_folder))
saveas(gcf, sprintf('%s/cs_magnitudeOrder', fig_folder), 'epsc')

%% 

[Mir_cs, Mi_cs] = MIR(Orbits);          % Mir of all data.
[Mir_cs2, Mi_cs2] = MIR(Orbits(1:6e4,:));  % Mir till the synchronization state.

for i=1:M
   Mir_cs(i,i) = NaN;
   Mir_cs2(i,i) = NaN;
end

Mir_cs = Mir_cs+ Mir_cs';
Mir_cs2 = Mir_cs2+ Mir_cs2';
%
figure
subplot(1,2,1)
h1 = imagesc(Mir_cs);
% c1.Limits = [0,0.2]; 
pbaspect([1 1 1])
%h1 = heatmap(Mir_cs, 'FontSize',20);
title('(a)', 'FontSize',30, 'Interpreter','latex')
% h1.NodeChildren(3).Title.Interpreter = 'latex';
ax= gca; 
ax.FontSize = 30; 

subplot(1,2,2)
h2 = imagesc(Mir_cs2);
pbaspect([1 1 1])
c1 = colorbar; 
ax= gca; 
ax.FontSize = 30; 
caxis([min([h1.CData(:); h2.CData(:)]), max([h1.CData(:); h2.CData(:)])]);

% h2 =heatmap(Mir_cs2, 'FontSize',20);
title('(b)', 'FontSize',30, 'Interpreter','latex')
% h2.NodeChildren(3).Title.Interpreter = 'latex';
% 
% set(gcf,'color','w'); % white background
% set(gcf, 'Position', [100 100 1600 600])

%%
saveas(gcf, sprintf('%s/MIR_cs.fig', fig_folder))
saveas(gcf, sprintf('%s/MIR_cs', fig_folder), 'epsc')
 
%% Rossler_phase Synchronisation

% Sync Measures for phases


Orbits = readmatrix('Rossler_Info/Orbits_rossler.txt');
tvalues = readmatrix('Rossler_Info/tvals_rossler.txt');
LYA = readmatrix('Rossler_Info/LYA.txt');
inst_phases = readmatrix('Rossler_Info/inst_phases.txt'); 
rossler_inf = readmatrix("Rossler_Info/info_0.020.txt");
discarded = rossler_inf(1,2);
coupling_strength = 0.02;
a = 0.15; f= 0.2; c= 10;
dw = 0.015; % frequency mismatch
w = linspace(1-dw, 1+dw,5); % 1+2*rand(M,1)*dw-dw;

t0 = 0;
tf = 1e4;
k=1e5;
dt = (tf-t0)/(k-1);
M = size(Orbits,2)/3;
%%


pws_mat = pairwise_sync(t0,tf,dt,inst_phases(:,1:5));
[corr_mat, pvals_dat] = corr(Orbits(:,1:5));
corr_mat = abs(corr_mat);

phasync_destroyed = phaseran(Orbits(:,1:5),1);
phasync_destroyed_pwsync = pairwise_sync(t0,tf-dt,dt,phasync_destroyed);
[phasyncdestroyed_signalcorr, pvals_twsd] = corr(phasync_destroyed);
phasyncdestroyed_signalcorr = abs(phasyncdestroyed_signalcorr);

for i=1:M
    pws_mat(i,i) = 0;
    corr_mat(i,i) = 0;
    pvals_twsd(i,i) =0;
    pvals_dat(i,i) = 0;
    phasync_destroyed_pwsync(i,i) = 0;
    phasyncdestroyed_signalcorr(i,i) = 0;
end
% Kuramoto phase order

Kuramoto_rho = global_sync(inst_phases);


%%
figure
subplot(3,3,1:3)
plot(tvalues(discarded+1:end), Kuramoto_rho, 'LineWidth',2)
xlabel('Time', 'interpreter', 'latex', 'Fontsize', 30)
ylabel('$\rho$', 'interpreter', 'latex', 'Fontsize', 30)
title('(a)', 'interpreter', 'latex', 'Fontsize', 30)
ax = gca;
ax.FontSize = 30; 

subplot(3,3,4)
h1 = imagesc(pws_mat);
title('(b)', 'FontSize',30, 'interpreter', 'latex')
ax = gca;
ax.FontSize = 30; 
pbaspect([1 1 1])
yticks([1,2,3,4,5])

subplot(3,3,7)
h4 =imagesc(phasync_destroyed_pwsync);
% title('', 'FontSize',30, 'interpreter', 'latex')
pbaspect([1 1 1])
ax = gca;
ax.FontSize = 30; 
colorbar;
caxis([min([h1.CData(:); h4.CData(:)]), max([h1.CData(:); h4.CData(:)])]);
xticks([1,2,3,4,5])
yticks([1,2,3,4,5])

%
subplot(3,3,5)
h2 = imagesc(corr_mat);
pbaspect([1 1 1])
title('(c)', 'FontSize',30, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30; 
yticks([1,2,3,4,5])


subplot(3,3,8)
h5 = imagesc(phasyncdestroyed_signalcorr);
pbaspect([1 1 1])
% title('', 'FontSize',30,'Interpreter','latex')
ax = gca;
ax.FontSize = 30; 
colorbar;
caxis([min([h2.CData(:); h5.CData(:)]), max([h2.CData(:); h5.CData(:)])]);
xticks([1,2,3,4,5])
yticks([1,2,3,4,5])


subplot(3,3,6)
h3 = imagesc(pvals_dat);
pbaspect([1 1 1])
title('(d)', 'FontSize',30, 'Interpreter','latex')
ax = gca;
ax.FontSize = 30; 
colorbar;
caxis([0,0.05]); %[min([h3.CData(:); h6.CData(:)]), max([h3.CData(:); h6.CData(:)])]);
yticks([1,2,3,4,5])

subplot(3,3,9)
h6 = imagesc(pvals_twsd);
pbaspect([1 1 1])
ax = gca;
ax.FontSize = 30; 
%title('', 'FontSize',30, 'Interpreter','latex')
colorbar;
caxis([0,0.05]); %[min([h3.CData(:); h6.CData(:)]), max([h3.CData(:); h6.CData(:)])]);
xticks([1,2,3,4,5])
yticks([1,2,3,4,5])



% set(gcf,'color','w'); % white background
% set(gcf, 'Position', [100 100 1600 1600])
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 16])  %% Size of the .png
%%
saveas(gcf, sprintf('%s/rossler_phase_sync_measures.fig', fig_folder))
saveas(gcf, sprintf('%s/rossler_phase_sync_measures', fig_folder), 'epsc')
%%
Dat = Orbits(:,1:5);
phasync_destroyed = phaseran(Orbits(:,1:5),1);


% Instantaneous phases and frequency from data and twinsd

inst_phases_dat = hilbert_phases(Dat);
inst_phases_twsd = hilbert_phases(phasync_destroyed);

inst_freq_dat = zeros(size(inst_phases_dat,1)-1, M);
inst_freq_twsd = zeros(size(inst_phases_twsd,1)-1, M);

for i=1:size(inst_phases_dat,1)-1
    if i > size(inst_phases_twsd,1)-1 
        inst_freq_dat(i,:) = (inst_phases_dat(i+1,:) -inst_phases_dat(i,:))/(dt*2*pi); 
    else
        inst_freq_dat(i,:) = (inst_phases_dat(i+1,:) -inst_phases_dat(i,:))/(dt*2*pi); 
        inst_freq_twsd(i,:) = (inst_phases_twsd(i+1,:) - inst_phases_twsd(i,:))/(dt*2*pi); 
    end
end

%%

[MIR_dat_orbit, ~] = MIR(Dat);
[MIR_twsd_orbit, ~] = MIR(phasync_destroyed);

[MIR_dat_phases, ~] = MIR(inst_phases_dat);
[MIR_twsd_phases, ~] = MIR(inst_phases_twsd);

[MIR_dat_freq, ~] = MIR(inst_freq_dat);
[MIR_twsd_freq, ~] = MIR(inst_freq_twsd);


for i=1:M
    MIR_dat_orbit(i,i) =0;
    MIR_twsd_orbit(i,i) = 0;
    MIR_dat_phases(i,i) =0;
    MIR_twsd_phases(i,i) = 0;
    MIR_dat_freq(i,i) = 0;
    MIR_twsd_freq(i,i) = 0;
end

MIR_dat_orbit = MIR_dat_orbit + MIR_dat_orbit';
MIR_twsd_orbit = MIR_twsd_orbit + MIR_twsd_orbit';
MIR_dat_phases = MIR_dat_phases + MIR_dat_phases';
MIR_twsd_phases = MIR_twsd_phases + MIR_twsd_phases';
MIR_dat_freq = MIR_dat_freq + MIR_dat_freq';
MIR_twsd_freq = MIR_twsd_freq + MIR_twsd_freq';
%%
figure
subplot(1,2,1)
h7 = imagesc(MIR_dat_orbit); 
title('(a)', 'FontSize',30, 'Interpreter','latex')
pbaspect([1 1 1])
ax= gca; 
ax.FontSize = 30;

subplot(1,2,2)
h8 = imagesc(MIR_twsd_orbit);
title('(b)','FontSize',30, 'Interpreter','latex')
colorbar;
caxis([min([h7.CData(:); h8.CData(:)]), max([h7.CData(:); h8.CData(:)])]);
pbaspect([1 1 1])
ax= gca; 
ax.FontSize = 30;
%%

saveas(gca, sprintf('%s/Rossler_twsd_MIR.fig', fig_folder))
saveas(gca, sprintf('%s/Rossler_twsd_MIR', fig_folder), 'epsc')

%%
% 
% subplot(3,2,3)
% h9 = imagesc(MIR_dat_phases);
% title('(c)','FontSize',30, 'Interpreter','latex')
% pbaspect([1 1 1])
% 
% subplot(3,2,4)
% h10 = imagesc(MIR_twsd_phases);
% title('(d)', 'FontSize',30, 'Interpreter','latex')
% colorbar;
% caxis([min([h9.CData(:); h10.CData(:)]), max([h9.CData(:); h10.CData(:)])]);
% pbaspect([1 1 1])
% 
% %%
% subplot(3,2,5)
% h10 = heatmap(MIR_dat_freq, 'Fontsize', 20);
% title('MIR for frequency')
% h10.NodeChildren(3).Title.Interpreter = 'latex';
% 
% subplot(3,2,6)
% h10 = heatmap(MIR_twsd_freq, 'Fontsize', 20);
% title('MIR for TWSD frequency')
% h10.NodeChildren(3).Title.Interpreter = 'latex';
% 
% 
% 
% %% best_case_TPR_FPR(similarity_matrix,A, plot_inf_net,plot_Roc_plot, barplot, label)
% 
% Adj  = [0 1 0 0 0;
%         1 0 1 1 0;
%         0 1 0 0 0;
%         0 1 0 0 1;
%         0 0 0 1 0];
% 
% M = size(Adj,2);
% 
% [TPR_FPR_ROC1, TPR_FPR_ROC_best1]  = best_case_TPR_FPR(MIR_dat_orbit,Adj, 0, 0, 0, 'Mir for Orbits');
% [TPR_FPR_ROC2, TPR_FPR_ROC_best2]  = best_case_TPR_FPR(MIR_twsd_orbit,Adj, 0, 0, 0, 'Mir for Orbits');
% [TPR_FPR_ROC3, TPR_FPR_ROC_best3]  = best_case_TPR_FPR(MIR_dat_phases,Adj, 0, 0, 0, 'Mir for Orbits');
% [TPR_FPR_ROC4, TPR_FPR_ROC_best4]  = best_case_TPR_FPR(MIR_twsd_phases,Adj, 0, 0, 0, 'Mir for Orbits');
% [TPR_FPR_ROC5, TPR_FPR_ROC_best5]  = best_case_TPR_FPR(MIR_dat_freq, Adj, 0, 0, 0, 'Mir for Orbits');
% [TPR_FPR_ROC6, TPR_FPR_ROC_best6]  = best_case_TPR_FPR(MIR_twsd_freq,Adj, 0, 0, 0, 'Mir for Orbits');
% FPR_list = zeros(6,1);
% TPR_list = zeros(6,1);
% 
% legend_list = {'Signal', 'TWSD Signal', 'Phase', 'TWSD Phases', 'Frequency', 'TWSD Frequency'};
% %%
% 
% figure
% plot(TPR_FPR_ROC_best1.FPR, TPR_FPR_ROC_best1.TPR , '*', 'MarkerSize', 20, 'Color', '#ffcc66');
% hold on
% plot(TPR_FPR_ROC_best2.FPR, TPR_FPR_ROC_best2.TPR, '*', 'MarkerSize', 20, 'Color', '#9900cc');
% alpha(.2)
% plot(TPR_FPR_ROC_best3.FPR, TPR_FPR_ROC_best3.TPR, '.', 'MarkerSize', 20, 'Color', '#ffcc66');    
% plot(TPR_FPR_ROC_best4.FPR, TPR_FPR_ROC_best4.TPR , '.', 'MarkerSize', 20, 'Color', '#9900cc');
% 
% plot(TPR_FPR_ROC_best5.FPR, TPR_FPR_ROC_best5.TPR, 'o', 'MarkerSize', 20, 'Color', '#ffcc66');
% plot(TPR_FPR_ROC_best6.FPR, TPR_FPR_ROC_best6.TPR, 'o', 'MarkerSize', 20, 'Color', '#9900cc');
% 
% plot(linspace(0,1,1e3), linspace(0,1,1e3), 'r--')
% xlabel('FPR', 'interpreter', 'latex', 'Fontsize', 20)
% ylabel('TPR', 'interpreter', 'latex', 'Fontsize', 20)
% title('ROC curve', 'interpreter', 'latex', 'Fontsize', 20)
% 
% legend(legend_list, 'interpreter', 'latex', 'Fontsize', 20, 'location', 'southeast')
% xlim([0 1])
% ylim([0 1])
% pbaspect([1 1 1])
% 
% %%
% saveas(gcf, sprintf('%s/ROC_orbits_twsd.fig', fig_folder))
% saveas(gcf, sprintf('%s/ROC_orbits_twsd.png', fig_folder))
% close


























