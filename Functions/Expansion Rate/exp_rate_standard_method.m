function [estimated_largest_expansion_rate, T] = exp_rate_standard_method(Dat, grid, T)


% The function estimate expansion rate from provided time series
% of 2 dimensional system using the idea of expansion rate. For higher 
% dimensional systems, make 'partition' suitable.
% (location function just taking the first 2 columns even if it has more.)
% but it would take longer time.

% cost-effective algorithm is needed to compute distance between pairs,
% look at the idea in Wolf-algorithm! (another partition to eliminate 
% the points located close to each other.)

niter = size(Dat,1);
[coordinates, ppb] = location(Dat,grid);

Nc = sum(ppb(:,2) >= niter/grid^2); % number of boxes have been occupied by at least 2.
numb_boxes = sort(ppb(ppb(:,2) >= niter/grid^2,1)); % idx of the occupied boxes.

%
% for i=1:unique(coordinates(:,end))
%     idx = coordinates(coordinates(:,end) == i, :);
%     
% 
% 
% end

if exist('T', 'var')
    T=T;
else
T = corr_decay(Dat,grid);
end

sum1 = 0;

for k = 1:Nc
        box = numb_boxes(k); 
        idx = coordinates(coordinates(:,end) == box, 1);
        pointsbox1 = Dat(idx, :);
        npts = size(pointsbox1,1);

        % delta = maxdist(pointsbox1);
        
        dist = zeros(npts*(npts-1)/2,1);
        counter = 0;
        for i=1:npts
            for j=i+1:npts
                counter = counter+1;
                dist(counter) = pdist2(pointsbox1(i,:), pointsbox1(j,:));
            end
        end
        
        delta = max(dist); % maximum distance between pairs in the initial cell.
        %
        idx = idx(idx+T < niter+1);
        iterated_points = Dat(idx+T, :);
        

       % DELTA = maxdist(iterated_points);
        
        npts2 = length(iterated_points);
        dist2 = zeros(npts2*(npts2-1)/2,1);
        counter = 0;
        for i=1:npts2
            for j=i+1:npts2
                counter = counter+1;
                dist2(counter) = pdist2(iterated_points(i,:), iterated_points(j,:));
            end
        end
        
        DELTA = max(dist2); %% maximum distance between pairs after iteration

        sum1 = sum1+log(DELTA/delta);

end

estimated_largest_expansion_rate = (1/T)*(sum1/Nc);

end
