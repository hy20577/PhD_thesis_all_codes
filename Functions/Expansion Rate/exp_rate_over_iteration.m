

function [exp_rate_T_over_grids, exp_rate] = exp_rate_over_iteration(Dat, expansion_err_rate)
% Dat here n by 2 matrix

%%%%%   Outputs
 
% exp_rate_T_over_grids
    % 1st column  --> Number of grids
    % 2nd column  --> Average Correlation Decay Time (T) from average of each cell iteration
    % 3rd column  --> Correlation Decay Time T from estimation from expansion rate (formula in paper.)
    % 4rd olumn  --> Correlation Decay Time T from itinerary network method
    % 5th column --> Expansion Rate from 2nd column.
    % 6th column --> exp rate from itinerary network!!


if ~exist('expansion_err_rate', 'var')
expansion_err_rate = 0.1;
end

%%%% 0-1 Scale the data  %%%%%

Dat(:,1) = (Dat(:,1) - min(Dat(:,1)))/(max(Dat(:,1) - min(Dat(:,1)))); 
Dat(:,2) = (Dat(:,2) - min(Dat(:,2)))/(max(Dat(:,2) - min(Dat(:,2)))); 

%%%%


niter = size(Dat,1);
[unbiased_grid_size_min, unbiased_grid_size_max] = Grid(Dat);
grid_range = unbiased_grid_size_min:unbiased_grid_size_max; % The grids satisfying the inequality 17 have been considered.
exp_rate_T_over_grids = zeros(length(grid_range),5);
sum2 = 0;

%setting up grid size randomly 

%idx_random_grid = randi(length(grid_range), 1);
rand_grid = 6; %grid_range(idx_random_grid);

for grid = grid_range(1):grid_range(end) % MAIN LOOP over grid size

        [coordinates, ppb] = location(Dat,grid);
        Nc = sum(ppb(:,2) >= niter/grid^2); % number of well-occupied cells, number of pts in the cell should be bigger than (or equal to) total number of pts/number of cells.
        numb_boxes = sort(ppb(ppb(:,2) >= niter/grid^2, 1)); % idx of the occupied boxes.


  
        sum1 = 0;
        sumT = 0;
        sumspecificied_cell_er= 0;
                    for k = 1:Nc % taking average over occupied cells.
                            box = numb_boxes(k); 
                            idx = coordinates(coordinates(:,end) == box, 1);
                            
                            pointsbox1 = Dat(idx, :);                    
   % estimation of delta
                            [delta, ~] = maxdist(pointsbox1);
                             
                            %if ~exist(number_of_iteration, 'var')
                                    DELTA = 0 ;
                                    T =0;

                                    while DELTA < sqrt(2)*(1-expansion_err_rate) 
                                    T= T+1;
                                    idx = idx(idx+T < niter+1);
                                    iterated_points = Dat(idx+T, :);% estimation of DELTA  
                                   [DELTA, ~] = maxdist(iterated_points);

%                                    if grid == rand_grid
%                                           exp_rate_on_specied_grid_over_time(T,1) = T;
%                                           exp_rate_on_specied_grid_over_time(T,2) = (1/T)*log(DELTA/delta);
%                                     end

                                    if T == 40
                                        break
                                    end

                                   end
                     sumT = sumT+T;
                            %else

                    sum1 = sum1+ (1/T)*log(DELTA/delta);
                   end
   
   exp_rate_T_over_grids(grid - grid_range(1)+1,1) = grid;
 
   %% Correlation Decay Time T from different estimation
   exp_rate_T_over_grids(grid - grid_range(1)+1, 2) = sumT/Nc; % Averaege T over cells
   exp_rate_T_over_grids(grid - grid_range(1)+1, 3) = (1/((Nc/sumT)*sum1/Nc))*log(grid); % T estimation from the e1
  
   exp_rate_T_over_grids (grid - grid_range(1)+1, 5) = (Nc/sumT)*sum1/Nc;  %% Exp rate from 
   

   sum2 = sum2+sum1/Nc;
   
end
    [estimated_largest_expansion_rate, ~, number_of_iter] = exp_rate_LYA(Dat);
    exp_rate_T_over_grids(:, 4) = number_of_iter(:,2);
    exp_rate =  sum2 / length(grid_range); 
    exp_rate_T_over_grids (:, 6) = estimated_largest_expansion_rate(:,2); % exp rate from itinerary network!!
   

%% Plot
% figure
% plot(exp_rate_on_specied_grid_over_time(:,1), exp_rate_on_specied_grid_over_time(:,2), '*--')
% xlabel('Number of Iterations', 'Interpreter','latex', 'FontSize',20)
% ylabel('Expansion Rate', 'Interpreter','latex', 'FontSize',20)
% title(sprintf('Expansion Rate Time Evolution on the %d by %d Partition', rand_grid,rand_grid),'Interpreter','latex', 'FontSize',20);
% yline(estimated_largest_expansion_rate(estimated_largest_expansion_rate(:,1)==rand_grid,2));
% legend('Expansion Rate Over Time', 'Estimated Average Exp Rate from itinerary networkß')

end

