%% Lyapunov Exponent for Ode systems

function [tvalues, Lyapunov, Orbits] = LyapunovforODEs_defined_length_deviation(x0,t0, tf, k, ext_odefun, len_of_deviation_vector)

%%% ode45 function is used to solve ode equation, needed to revise for
%%% other solution.

% ext_odefun - extended ode system needed as input
% x0: initial condition for the ode system
% t0: initial time
% tf: final time
% k: number of iteration
% How many data you want to discard in computation of LE as 

n = length(x0); % length of original ode
stept = (tf-t0)/k;
sum1= 0;

% Setting initial conditions
t= t0;
y0 = zeros(1, n+len_of_deviation_vector);
y0(1:n) = x0;

M = eye(sqrt(len_of_deviation_vector));

y0(n+1:end) = reshape(M, 1, len_of_deviation_vector);

Orbits = zeros(k,n);
tvalues = zeros(k,1);

    for ITERLAP=1:k   %Main loop to compute Lyapunov
        
        tspan= [t t+stept];
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) ext_odefun(t,y), tspan, y0, options);
        tvalues(ITERLAP) = t+stept;
        Y = Y(end,:);
        y0(1:n) = Y(1:n);
        Orbits(ITERLAP,:) = Y(1:n);
        
        K = reshape(Y(n+1:end),sqrt(len_of_deviation_vector),sqrt(len_of_deviation_vector));
        
%         for j=1:n
%      
%             K(:,j) = Y(n*j+1 : n*(j+1));  %% create deviation matrix
%             
%         end
        
        K = gramschmidt(K);    %% Orthogonalization
        
        norm_K = zeros(1, sqrt(len_of_deviation_vector));
        NK = zeros(sqrt(len_of_deviation_vector), sqrt(len_of_deviation_vector));
        
        for j=1:sqrt(len_of_deviation_vector)
            norm_K(j) = norm(K(:,j));
            NK(:,j) = K(:,j)/norm_K(j);
        end
        
        sum1 = sum1 + log(norm_K);
       
        t = t+stept;
        
        y0(n+1:end) = reshape(NK,1,(len_of_deviation_vector));
        
%         for j=1:n
%          y0(n*j+1 : n*(j+1)) = K(:,j)/norm_K(j);  %% prepare deviation matrix for next iteration.   
%         end
   
       if mod(ITERLAP,100) == 0
            fprintf('Progression: %s %% \n', num2str(ITERLAP*100/k))
       end
   end
   Lyapunov = sum1/(tf-t0);
   Lyapunov = sort(Lyapunov,'descend'); 

end

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
        c = c - (a*a'*b)/(a'*a);
    end
    U(:,ii) = c;
end
end

