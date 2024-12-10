function Adj = ringAdjMat(M,r)

if ~isinteger(r)
    warning(sprintf('Number of neighbour is not integer. Assumed r = %d instead r=%.2f', ceil(r), r))
end

r = ceil(r);

    Adj = zeros(M,M);    
    for i=1:M
        for j=1:M
            if abs(mod(i,M)-mod(j,M)) < r+1 && i~=j
                Adj(i,j)=1;
            elseif abs(mod(i,M)-mod(j,M)) > M-r-1 && i~=j
                Adj(i,j)=1;
            end      
      
        end
    end
end