function [Orbits, Lya_overtime] = Circle_generator(Adjacency, alpha, x0, r, niter, time_window, std_threshold, max_discard, savefigs, label)

M = size(Adjacency,2);
if length(r) ==1
    r = ones(M,1)*r;
end


if max_discard < time_window
    max_discard = time_window;
    warning('max_discard = time_window set up automatically')
end


k = length(x0);
w = eye(k);
v = x0;
sum1 = 0;
Orbits = zeros(niter, k);
    for i = 1:time_window
        K = zeros(k,k);   % matrix stores deviation vectors in columns.
            J = Jacob_cm(v, Adjacency, alpha); 

            for j =1:k
            K(:,j) = J*w(:,j);  % images of w() under Jacobian matrix.  
            end
            
            ort_K = gramschmidt(K);  %% gs function produces orthonormal vectors.
            y = zeros(1,k);
            
            for j=1:k
            y(j) = norm(ort_K(:,j));
                if sum(ort_K(:,j)) ~= 0 
                    w(:,j) = ort_K(:,j)/norm(ort_K(:,j));
                else
                    w(:,j) = 0;
                end
            end
            sum1 = sum1 + log(y);
            Lya_overtime(i,:) = sort(sum1/i, 'descend');
            v= coupled_circle_fnc(v,r,Adjacency,alpha);
    end

     sum2 = 0;

    while std(Lya_overtime(end-time_window+1:end,1)) > std_threshold 
        sum2 = sum2+1;
        K = zeros(k,k);   % matrix stores deviation vectors in columns.
        J = Jacob_cm(v, Adjacency, alpha); 
                for j =1:k
                K(:,j) = J*w(:,j);  % images of w() under Jacobian matrix.  
                end
            ort_K = gramschmidt(K);  %% gs function produces orthonormal vectors.
            y = zeros(1,k);

            for j=1:k
            y(j) = norm(ort_K(:,j));
                if sum(ort_K(:,j)) ~= 0 
                    w(:,j) = ort_K(:,j)/norm(ort_K(:,j));
                else
                    w(:,j) = 0;
                end
            end

            sum1 = sum1 + log(y);
            Lya_overtime(time_window+sum2,:) = sort(sum1/(sum2+time_window), 'descend');
            v= coupled_circle_fnc(v,r,Adjacency,alpha);

            if sum2 > max_discard-time_window 
                warning(sprintf('std_threshold is too small to stabilize; first %d discarded', time_window+max_discard))
                break
            end

    end
    
    for i = 1:niter
   
            Orbits(i,:) = v;
            K = zeros(k,k);   % matrix stores deviation vectors in columns.
            J = Jacob_cm(v, Adjacency, alpha); 

            for j =1:k
            K(:,j) = J*w(:,j);  % images of w() under Jacobian matrix.  
            end
            
            ort_K = gramschmidt(K);  %% gs function produces orthonormal vectors.
            y = zeros(1,k);
            
            for j=1:k
            y(j) = norm(ort_K(:,j));
                if sum(ort_K(:,j)) ~= 0 
                    w(:,j) = ort_K(:,j)/norm(ort_K(:,j));
                else
                    w(:,j) = 0;
                end
            end
  
            sum1 = sum1 + log(y);
           
            Lya_overtime(i+time_window+sum2,:) = sort(sum1/(i+time_window+sum2), 'descend'); 
            
            v= coupled_circle_fnc(v,r,Adjacency,alpha);
            
        if  mod(i,1e3) == 0
            fprintf('Completed job : %s %% \n', floor(num2str(i*100/niter)))
        end
        
        
    end
if savefigs == 1
    currpath = pwd; 
    
    if ~exist(sprintf('%s/Circle_Info_%s', currpath, label), 'dir')
        mkdir(sprintf('%s/Circle_Info_%s', currpath, label))
    else
        rmdir(sprintf('%s/Circle_Info_%s', currpath, label), 's')
        mkdir(sprintf('%s/Circle_Info_%s', currpath, label))
    end
    
    folder_path = sprintf('%s/Circle_Info_%s', currpath, label); 

    info = {'discarded', time_window+sum2; 'MLE', Lya_overtime(end,1); 'niter', niter;'std_threshold', std_threshold; ...
      'coupling strengtht', alpha; 'time window', time_window; 'max_discard', max_discard};
    writecell(info, sprintf('%s/info_%.2f.txt',folder_path, alpha))
    
    writematrix(Lya_overtime(:,1), sprintf('%s/MLE_%.2f.txt',folder_path, alpha))
    
    discard = length(Lya_overtime(:,1)) - niter;
    figure;
    loglog(Lya_overtime(:,1))
    xlim([0 length(Lya_overtime(:,1))])
    ylim([-abs(min(Lya_overtime(:,1))), 1.5*max(Lya_overtime(:,1))])
    xline(discard, 'r--')
    yline(0, 'r-')
    hold on
    patch([0 discard discard 0], [-abs(min(Lya_overtime(:,1))) -abs(min(Lya_overtime(:,1))) 1.5*max(Lya_overtime(:,1)) 1.5*max(Lya_overtime(:,1))],'red', 'FaceAlpha', 0.3)
    t = text((discard)/2,-abs(min(Lya_overtime(:,1)))/2,'Transient Period');
    t.Rotation = 90;
    hold off
    saveas(gca, sprintf('%s/MLE_over_Time_%.2f.png',folder_path, alpha))
    close all
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%


function U = gramschmidt(V)
% GRAM_SCHMIDT - Classic Gram-Schmidt Process
%   Input vectors should be the columns of input matrix.
%   Output = unitary, orthogonal vectors in columns of matrix

n = size(V, 2); % number of columns

U(:,1) = V(:,1);

% find next orthogonal column
    for ii = 2:n
        b = V(:,ii);
        c = b;
        for k = 1:ii-1
            a = U(:,k);
            if a ~= 0 
                c = c - (a*a'*b)/(a'*a);
            end
        end
        U(:,ii) = c;
    end
end

function Circle = coupled_circle_fnc(x,r ,A, alpha)
% coupled_circular_fnc gives the iterations of circular map
%%% To compute the Lyapunov Exponents with the LyapunovforMAPS 
%%% functions, it requires to pass variables except 'x' as a global.  
M= length(A(:,1));   %% Number of nodes, in other words, length of A.
if length(r) == 1
    r = ones(M,1)*r;
end

Circle = zeros(M,1);    
    for i=1:M %rows loop
        sum2=0;
        for j=1:M %columns loop
            sum2 = sum2 + A(i,j) * circ(r(j), x(j));
        end
        if sum2 ~= 0
            coupling_term = sum2 * (alpha/sum(A(i,:)));
        else
            coupling_term = 0;
        end
            Circle(i) = circ(r(i), x(i)) * (1-alpha) + coupling_term;
    end
    
    function circ = circ(r, x)
       K= 6.9115; 
       circ = x + r - (K/(2*pi))* sin(2*pi*x);
       circ = mod(circ,1);
    end

end

function Jacob_cm= Jacob_cm(x, A, alpha)
%Jacob_cm calculates the jacobian matrix in the given set of points
% x is vector of initial points.

k = length(x);
Jacob_cm = zeros(k,k);
    for i=1:k
        for j=1:k
            if i == j
                Jacob_cm(i,j) = df_cm(x(i))*(1-alpha);
            else
                if sum(A(i,:) ~= 0)
                    Jacob_cm(i,j) = (alpha/sum(A(i,:)))*A(i,j)*df_cm(x(j));
                else
                    Jacob_cm(i,j) = 0;  
                end
            end
        end
    end
    
    function df_cm = df_cm(x)
    K = 6.9115;
    df_cm = 1 - K* cos(2*pi*x);
    end
end

