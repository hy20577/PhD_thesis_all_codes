function [TPR_mat_heatmap, FPR_mat_heatmap, inf_nets, Mir_overall, output_all] = heatmap_new_method_general(Orbits, p_sl, Adj, varargin)

if nargin == 3
    pc1 = 0:0.1:1;
    pc2 = 0:0.1:1;
    label = '';
elseif nargin == 4
    pc1 =  varargin{1};
    pc2 = 0:0.1:1;
    label ='';
elseif nargin == 5
    pc1 = varargin{1};
    pc2 = varargin{2};
    label ='';
elseif nargin == 6
    pc1 = varargin{1};
    pc2 = varargin{2};
    label =varargin{3};
end

time_range = size(Orbits,1);
M = size(Adj,2);

l1 = length(pc1);
l2 = length(pc2);
output_all = zeros(l1*l2, 5); % pc1, pc2, TPR, FPR, distance  
TPR_matrix = zeros(l1,l2);
FPR_matrix = zeros(l1,l2);
dens = zeros(l1,l2);
dens_orig = sum(sum(Adj))/(M*(M-1));
counter =0;

%[Mir, ~] = MIR(Orbits);
% sumdiff = zeros(l1*l2,3);
inf_nets = zeros(M,M,l1*l2);
% Mir = Mir + Mir';
Mir_overall = zeros(1/p_sl+1,M*(M-1)/2,l1*l2);
for i=1:l1
    for j=1:l2
        %Laplacian_mat= zeros(M,M);
        counter = counter+1;
        [MIR_all, adj_MIR_fdr, TP_FP_Rate] = Mir_surrogate_data(Orbits, p_sl, time_range, "HT_mix_A_theta", Adj,'', pc1(i), pc2(j));
        Mir_overall(:,:,counter) = MIR_all;
        output_all(counter,1) = pc1(i);
        output_all(counter,2) = pc2(j);
        output_all(counter,3:4) = TP_FP_Rate(2:3);
        output_all(counter,5) = sqrt((1-TP_FP_Rate(2))^2+TP_FP_Rate(3)^2);

        TPR_matrix(i,j) = TP_FP_Rate(2);
        FPR_matrix(i,j) = TP_FP_Rate(3);
        dens(i,j) = sum(sum(adj_MIR_fdr))/(M*(M-1)) - dens_orig; % density of the inferred network - original network.
        
        inf_nets(:,:,counter) = adj_MIR_fdr;

    % Computing MIR values for original timeseries and Laplacian-generated
     % if all(sum(adj_MIR_fdr,1) ~= 0)  %% Cannot include the isolated nodes !!! MIR function doesnot work for isolated one, 
       
        sprintf('%d out of %d has done', counter, l1*l2)
    end
end
       output_all = sortrows(output_all,5);

       %% Create directory
% dir_name = sprintf('%s/%s_heatmaps_%.3f', pwd,label, p_sl);
% if ~exist(dir_name, 'dir')
%     mkdir(dir_name)
% else
%     rmdir(dir_name,'s')
%     mkdir(dir_name)
% end

% %%Change matrix structure for heatmap

TPR_mat_heatmap = zeros(l1,l2);
FPR_mat_heatmap = zeros(l1,l2);
% dens_mat_heatmap = zeros(l1,l2);

for i=1:l1
TPR_mat_heatmap(i,:) = TPR_matrix(l1-i+1, :);
FPR_mat_heatmap(i,:) = FPR_matrix(l1-i+1, :);
% dens_mat_heatmap(i,:) = dens(l1-i+1, :);
end


%   %% Heatmaps
% figure
% subplot(1,3,1)
% h1 = heatmap(TPR_mat_heatmap, 'Colormap',parula);
% title('TPR')
% xlabel('Strength of Phase Randomisation');
% ylabel('Strength of Amplitude Randomisation')
% h1.XData = string(pc2);
% h1.YData = string(linspace(pc1(end),pc1(1),l1));
% h1.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h1.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h1.NodeChildren(3).Title.Interpreter = 'latex';
% h1.ColorLimits = [0 1]; 
% 
% subplot(1,3,2)
% h2 = heatmap(FPR_mat_heatmap, 'Colormap',parula);
% ylabel('Strength of Amplitude Randomisation');
% xlabel('Strength of Phase Randomisation')
% title('FPR')
% h2.XData = string(pc2);
% h2.YData = string(linspace(pc1(end),pc1(1),l1));
% h2.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h2.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h2.NodeChildren(3).Title.Interpreter = 'latex';
% h2.ColorLimits = [0 1];
% 
% subplot(1,3,3)
% h3 = heatmap(dens_mat_heatmap, 'Colormap', parula);
% title('Density of Inferred Network')
% ylabel('Strength of Amplitude Randomisation');
% xlabel('Strength of Phase Randomisation')
% h3.XData = string(pc2);
% h3.YData = string(linspace(pc1(end),pc1(1),l1));
% 
% h3.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h3.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h3.NodeChildren(3).Title.Interpreter = 'latex';

%% Figures
% figure
% compactheatmap_fromTFPR(TPR_mat_heatmap,FPR_mat_heatmap, label)
% saveas(gcf, sprintf('%s/Heatmaps_TPR_FPR.fig',dir_name))
% saveas(gcf, sprintf('%s/Heatmaps_TPR_FPR.pdf', dir_name))
% 
% writematrix(TPR_mat_heatmap, sprintf('%s/TPR_mat.txt', dir_name))
% writematrix(FPR_mat_heatmap, sprintf('%s/FPR_mat.txt', dir_name))
% writematrix(dens_mat_heatmap, sprintf('%s/dens_mat.txt', dir_name))
% writematrix(Mir_overall,sprintf('%s/Mir_overall.txt', dir_name) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ldat_MIRdiff = sumdiff(sumdiff(:,3)~=0,:);
% if type_similarity == "MIRdiff"
% proposed_pairs = Ldat_MIRdiff(Ldat_MIRdiff(:,3)== min(Ldat_MIRdiff(:,3)),:);
% elseif type_similarity == "MIRcorr"
% proposed_pairs = Ldat_MIRdiff(Ldat_MIRdiff(:,3)== max(Ldat_MIRdiff(:,3)),:);
% elseif type_similarity == "eigdiff"
% proposed_pairs = Ldat_MIRdiff(Ldat_MIRdiff(:,3)== min(Ldat_MIRdiff(:,3)),:);
% end
% %sumdiff = sortrows(sumdiff,3);
% 
% sprintf('Proposed network produced from pc1 = %d and pc2 = %d', proposed_pairs(1,1), proposed_pairs(1,2))
% 

%% The data at a first glance

% [Mir, ~] = MIR(Orbits);
% phases = angle(hilbert(Orbits));
% [Mir_amplitudes ,~]= MIR(abs(hilbert(Orbits)));
% pws = pairwise_sync_of_phases(phases);
% Mir = Mir + Mir';
% Mir_amplitudes = Mir_amplitudes+ Mir_amplitudes';
% 
% figure
% subplot(2,2,1)
% h4 = heatmap(abs(corr(Orbits)));
% title('PC of trajectory');
% h4.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h4.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h4.NodeChildren(3).Title.Interpreter = 'latex';
% 
% subplot(2,2,2)
% h5 = heatmap(Mir);
% title('Mir of Trajectory');
% h5.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h5.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h5.NodeChildren(3).Title.Interpreter = 'latex';
% 
% subplot(2,2,3)
% h6 = heatmap(Mir_amplitudes);
% title('Mir of Inst. Amplitudes');
% h6.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h6.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h6.NodeChildren(3).Title.Interpreter = 'latex';
% 
% 
% subplot(2,2,4)
% h7 = heatmap(pws);
% title('Pairwise Phase Synchronisation');
% h7.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
% h7.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
% h7.NodeChildren(3).Title.Interpreter = 'latex';
% 
% saveas(gcf, sprintf('%s/FirstLookatData.fig',dir_name))
% saveas(gcf, sprintf('%s/FirstLookatData.pdf', dir_name))

% writematrix(sumdiff, sprintf('%s/%s.txt', dir_name, label))

%% Important txt outputs
% writematrix(inf_nets, sprintf('%s/InferredNetworks.txt', dir_name))
% writematrix(output_all, sprintf('%s/output_all.txt', dir_name))

% writematrix(proposed_pairs, sprintf('%s/Proposed_Pairs.txt', dir_name))
close all;
end