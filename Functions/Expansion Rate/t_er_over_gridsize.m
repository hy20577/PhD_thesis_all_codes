function [t_N, av_exp_rate, T_exp_itinerary_network ] = t_er_over_gridsize(Dat, grid_range, expansion_err_rate)

% T correlation decay time estimation from expansion rate and itinerary
% network model!

% Dat: 2 columns timeseries. if it is more code only consider first 2.

% grid range: optional. If undefined, the range of grid size, in eq.17
% successfull network inference paper of Chris, will be considered.

% expansion_err_rate -> what rate of the plane should data spread over? 0.8
% or 0.9 common values.


%%%% 0-1 Scale the data  %%%%%

Dat(:,1) = (Dat(:,1) - min(Dat(:,1)))/(max(Dat(:,1) - min(Dat(:,1)))); 
Dat(:,2) = (Dat(:,2) - min(Dat(:,2)))/(max(Dat(:,2) - min(Dat(:,2)))); 


if ~exist('grid_range', 'var')
[unbiased_grid_size_min, unbiased_grid_size_max] = Grid(Dat);
grid_range = unbiased_grid_size_min:unbiased_grid_size_max; % The grids satisfying the inequality 17 have been considered.
end

if ~exist('expansion_err_rate', 'var')
expansion_err_rate = 0.1;
end

niter = size(Dat,1);
t_N = zeros(length(grid_range), 1);
counter = 0;
av_exp_rate = zeros(length(grid_range), 2);
T_exp_itinerary_network = zeros(length(grid_range),3);
CPU_time = zeros(length(grid_range), 1);
% t_N_expformula=  0;

for grid = grid_range(1):grid_range(end) % MAIN LOOP over grid size
start_time = tic(); 

        [coordinates, ppb] = location(Dat,grid);
        Nc = length(unique(coordinates(:,4)));  %% The number of occupied cell
        NcTilda = sum(ppb(:,2) >= niter/Nc); % Well occupied cells have been chosen if the number of points inside the cell is more than the mean number of pts over occupied cells.
        numb_boxes = sort(ppb(ppb(:,2) >= niter/Nc, 1)); % idx of the occupied boxes.
        sumT = 0;
        exp_rate_over_cells = 0; 
        exp_rate_itnet_over_cells = 0;
        T_itnet = corr_decay(Dat, grid); %% Correlation decay time from itinerary network method!
        
       parfor k = 1:NcTilda % taking average over occupied cells.   % parfor
                            box = numb_boxes(k); 
                            idx = coordinates(coordinates(:,end) == box, 1);
                            
                            pointsbox1 = Dat(idx, :);                    
                          
                            [delta, ~] = convhull_distance(pointsbox1);
                            %[delta, ~] = maxdist(pointsbox1);   % estimation of delta

                            
                            %% Also, computation of expansion rate from itinerary network
                                    idx_itnet = idx(idx+T_itnet < niter+1);
                                    iterated_points = Dat(idx_itnet+T_itnet, :);
                                   [DELTA_itnet, ~] = convhull_distance(iterated_points); % estimation of DELTA  


        exp_rate_itnet_over_cells = exp_rate_itnet_over_cells + (1/T_itnet)*log(DELTA_itnet/delta);

                            %% Iteration
     
                                    DELTA = 0 ;
                                    T =0;

                                    while DELTA < sqrt(2)*(1-expansion_err_rate) 
                                    T= T+1;
                                    idx = idx(idx+T < niter+1);
                                    iterated_points = Dat(idx+T, :);
                                   [DELTA, ~] = convhull_distance(iterated_points); % estimation of DELTA  

                                       if T == 100
                                            break
                                       end
                                    end                
            exp_rate_over_cells = exp_rate_over_cells + (1/T)*log(DELTA/delta);                        
            sumT = sumT + T;
           
       end  
    % t_N_expformula = t_N_expformula + (1/exp_rate_over_cells)*log(grid); 
     counter = counter +1 ;
     t_N(counter,1) = grid; 
     t_N(counter,2) = sumT/NcTilda;
      
     av_exp_rate(counter,1) = grid; 
     av_exp_rate(counter,2) = exp_rate_over_cells/NcTilda;


    T_exp_itinerary_network(counter, 1) = grid;
    T_exp_itinerary_network(counter, 2) =  T_itnet; 
    T_exp_itinerary_network(counter, 3) = exp_rate_itnet_over_cells/NcTilda;
    
    sprintf('Grid %d has done!', grid)
    CPU_time(counter) = toc(start_time);
    writematrix(CPU_time, 'CPU.txt')
end
    % t_N_expformula = t_N_expformula/length(grid_range);
end