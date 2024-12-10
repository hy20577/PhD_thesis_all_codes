

function [estimated_largest_expansion_rate, e1, number_of_iter] = exp_rate_LYA(Dat, grid_range,fig_switch)

%%%%% INPUTS 



%%% OPTIONAL

% grid_range is a length 2 vector describing min and max grid size. If not
    % defined, just considering the unbiased grid sizes.  
% 'fig_switch = 0' to not produce figures, otherwise you will get figures.

%%%% OUTPUTS



% The function estimate expansion rate from provided time series
% of 2 dimensional system using the idea of expansion rate. For higher 
% dimensional systems, make 'partition' suitable.
% (location function just taking the first 2 columns even if it has more.)
% but it would take longer time.

% cost-effective algorithm is needed to compute distance between pairs,
% look at the idea in Wolf-algorithm! (another partition to eliminate 
% the points located close to each other.)

%%%% 0-1 Scale the data  %%%%%

Dat(:,1) = (Dat(:,1) - min(Dat(:,1)))/(max(Dat(:,1) - min(Dat(:,1)))); 
Dat(:,2) = (Dat(:,2) - min(Dat(:,2)))/(max(Dat(:,2) - min(Dat(:,2)))); 
%%


niter = size(Dat,1);

counter = 0;
counter2 = 0;
[unbiased_grid_size_min, unbiased_grid_size_max] = Grid(Dat);

if ~exist('grid_range', 'var') % || grid_range == ""
    grid_range = unbiased_grid_size_min:unbiased_grid_size_max;
end

if ~exist("fig_switch", 'var')
    fig_switch = 0;
end


gno = randi([grid_range(1), grid_range(end)], 1);
e1 = zeros(length(grid_range),2);
estimated_largest_expansion_rate = zeros(length(unbiased_grid_size_min:unbiased_grid_size_max),2);

number_of_iter = zeros(length(grid_range),2);

for grid = grid_range(1):grid_range(end) % MAIN LOOP over grid size
        
        T = corr_decay(Dat,grid);   %%% Correlation Decay time T or T+1.

        counter = counter +1;

        number_of_iter(counter,1) = grid;
        number_of_iter(counter,2) = T;

        [coordinates, ppb] = location(Dat,grid);
        Nc = sum(ppb(:,2) >= 1); % number of occupied cells.
        Nctilda = sum(ppb(:,2) >= niter/Nc); % number of well-occupied cells. The cells which are occupied by at least 'niter/Nc' pts will be considered.
        numb_boxes = sort(ppb(ppb(:,2) >= niter/Nc, 1)); % idx of the welll-occupied boxes.
        idx_bno = randi([1, length(numb_boxes)],1);
        bno = numb_boxes(idx_bno);

        sum1 = 0;
                    for k = 1:Nctilda % taking average over occupied cells.
                            box = numb_boxes(k); 
                            idx = coordinates(coordinates(:,end) == box, 1);
                            
                            pointsbox1 = Dat(idx, :);                    
   % estimation of delta
                            [delta, nodes_small_delta] = maxdist(pointsbox1);

   % iterated points
                            idx = idx(idx+T < niter+1);
                            iterated_points = Dat(idx+T, :);
   % estimation of DELTA  
                           [DELTA, nodes_delta_capital] = maxdist(iterated_points);
                           
                             
                            if gno == grid && bno == box
                                  pts_coor = pointsbox1;
                                  iter_coor = iterated_points;
                                  coor_deltacapital = [box, nodes_delta_capital];  
                                  coor_deltasmall = [box, nodes_small_delta];
                            end

                           
                            sum1 = sum1+log(DELTA/delta);
                    
                    end
                    
                    e1(counter,1) =  grid;
                    e1(counter,2) = (1/T)*(sum1/Nc);


                 
                    if (grid > unbiased_grid_size_min || grid == unbiased_grid_size_min) && (grid < unbiased_grid_size_max || grid == unbiased_grid_size_max)
                        counter2 = counter2 + 1;
                        estimated_largest_expansion_rate(counter2,1) = grid;
                        estimated_largest_expansion_rate(counter2,2) = (1/T)*(sum1/Nc);
                    end

end

 %%%%%%    Plots   %%%%%%%%%%%

if fig_switch ~= 0
    
    figure
    plot(e1(:,1), e1(:,2), '*--')
    hold on
    plot(estimated_largest_expansion_rate(:,1), estimated_largest_expansion_rate(:,2), 'ro-')
    hold off
    xlabel("Grid Size", "Interpreter","latex", "FontSize",20)
    ylabel('Expansion Rate', "Interpreter","latex", "FontSize",20)
    title("Expansion Rate over grid size", "Interpreter","latex", "FontSize",20)
    
    figure
    plot(Dat(:,1), Dat(:,2), '.')
    xlabel("$x_{1}$", 'Interpreter','latex', 'FontSize', 20)
    ylabel("$x_{2}$", 'Interpreter','latex', 'FontSize', 20)
    title("Phase Space", 'Interpreter','latex', 'FontSize', 20)

    figure
    tiledlayout(1,2)
    nexttile
    plot(pts_coor(:,1), pts_coor(:,2), '.')
    xlim([0 1])
    ylim([0 1])
    hold on
    plot([coor_deltasmall(2), coor_deltasmall(4)], [coor_deltasmall(3), coor_deltasmall(5)], 'ro-', 'LineWidth',2)
    xlabel("$x_{1}$", 'Interpreter','latex', 'FontSize', 20)
    ylabel("$x_{2}$", 'Interpreter','latex', 'FontSize', 20)
    title("Initial cell - $\delta$ ", 'Interpreter','latex', 'FontSize', 20)
    for i=1:gno
        hold on
        xline(i*1/gno)
        yline(i*1/gno)
    end
    pbaspect([1 1 1])
    nexttile
    plot(iter_coor(:,1), iter_coor(:,2), '.')
    hold on
    plot([coor_deltacapital(2), coor_deltacapital(4)], [coor_deltacapital(3), coor_deltacapital(5)], 'ro-', 'LineWidth',2)
    xlabel("$x_{1}$", 'Interpreter','latex', 'FontSize', 20)
    ylabel("$x_{2}$", 'Interpreter','latex', 'FontSize', 20)
    title(sprintf("After %d Iteration - %s", number_of_iter(number_of_iter(:,1) == gno,2), '$\Delta$'),'Interpreter','latex', 'FontSize', 20)
    for i=1:gno
        hold on
        xline(i*1/gno)
        yline(i*1/gno)
    end
    pbaspect([1 1 1])

end


end