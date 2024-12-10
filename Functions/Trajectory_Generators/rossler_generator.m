function [tvalues, LYA, Orbits] = rossler_generator( x0, t0, tf, k, coupling_strength, a, f, c, w, Adj, time_window, std_threshold, max_discard)

%%% ode45 function is used to solve ode equation, needed to revise for
%%% other solution.

% ext_odefun - extended ode system needed as input
% x0: initial condition for the ode system
% t0: initial time
% tf: final time
% k: number of iteration
% a, f, and c model parameters; w internal frequency of the oscillators.
% Adj = Adjacecny matrix



if max_discard < time_window
    max_discard = time_window;
    warning('max_discard = time_window set up automatically')
end


n = length(x0); % dimension of extended ode function.
stept = (tf-t0)/k;
sum1= 0;

% Setting initial conditions
t= t0;
y0 = zeros(1, n^2+n);
y0(1:n) = x0;

for i=1:n
    y0((n+1)*i) = 1; % initial deviation vectors: unit vectors.
end

Orbits = zeros(k,n);
% tvalues = zeros(k,1);
% evolution_LYA = zeros(k,n+1);

for i = 1:time_window
    
        tspan= [t t+stept];
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) rossler_ext(t, y, Adj,w, a, f, c, coupling_strength), tspan, y0, options);

        tvalues(i) = t+stept;
        Y = Y(end,:);
        y0(1:n) = Y(1:n);

        K = reshape(Y(n+1:n^2+n),n,n);    
        K = gramschmidt(K);    %% Orthogonalization
        
        norm_K = zeros(1, n);
        NK = zeros(n,n);
        
        for j=1:n
            norm_K(j) = norm(K(:,j));
            NK(:,j) = K(:,j)/norm_K(j);
        end
       
        sum1 = sum1 + log(norm_K);
        t = t+stept;
        
        LYA(i,:) = sort(sum1/(t-t0), "descend");

        
        y0(n+1:n^2+n) = reshape(NK,1,n^2);
end

sum2 =0;
 while std(LYA(end-time_window+1:end,1)) > std_threshold 
        sum2 = sum2+1;
       
        tspan= [t t+stept];
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) rossler_ext(t, y, Adj,w, a, f,c, coupling_strength), tspan, y0, options);
        tvalues(i+sum2) = t+stept;
        Y = Y(end,:);
        y0(1:n) = Y(1:n);

        K = reshape(Y(n+1:n^2+n),n,n);    
        K = gramschmidt(K);    %% Orthogonalization
        
        norm_K = zeros(1, n);
        NK = zeros(n,n);
        
        for j=1:n
            norm_K(j) = norm(K(:,j));
            NK(:,j) = K(:,j)/norm_K(j);
        end

        sum1 = sum1 + log(norm_K);
        LYA(i+sum2,:)= sort(sum1/(t-t0), "descend");

        t = t+stept;
        y0(n+1:n^2+n) = reshape(NK,1,n^2);
     
    if mod(sum2,100) == 0
        fprintf('First %d data points have been removed. \n std = %.4f \n', sum2, std(LYA(end-time_window+1:end,1)))
    end

     if sum2 > max_discard-time_window 
        warning(sprintf('std_threshold is too small to stabilize; first %d discarded', max_discard))
        break
     end
 end

    for ITERLAP=1:k   %Main loop to compute Lyapunov
        
        tspan= [t t+stept];
        options = odeset('RelTol',1e-5);
        [T, Y] = ode45(@(t,y) rossler_ext(t, y, Adj,w, a, f,c, coupling_strength), tspan, y0, options);
        tvalues(time_window+sum2+ITERLAP) = t+stept;
        Y = Y(end,:);
        y0(1:n) = Y(1:n);
        Orbits(ITERLAP,:) = Y(1:n);
        
        K = reshape(Y(n+1:n^2+n),n,n);
        K = gramschmidt(K);    %% Orthogonalization
        
        norm_K = zeros(1, n);
        NK = zeros(n,n);
        
        for j=1:n
            norm_K(j) = norm(K(:,j));
            NK(:,j) = K(:,j)/norm_K(j);
        end
       
        sum1 = sum1 + log(norm_K);
        LYA(i+sum2+ITERLAP,:)= sort(sum1/(t-t0), "descend");

       % evolution_LYA(ITERLAP,1) = t; 
       % evolution_LYA(ITERLAP,2:end) = sum1/(t-t0);

        t = t+stept;
        
        y0(n+1:n^2+n) = reshape(NK,1,n^2);
        
%         for j=1:n
%          y0(n*j+1 : n*(j+1)) = K(:,j)/norm_K(j);  %% prepare deviation matrix for next iteration.   
%         end
   
    if mod(ITERLAP,1e2) == 0
        fprintf('Progression: %s %% \n', num2str(ITERLAP*100/k))
    end
    
    end
% Lyapunov_spectrum = sum1/(tf-t0);
% Lyapunov_spectrum = sort(Lyapunov_spectrum, 'descend');
% writematrix(evolution_LYA, 'Time_evolution_Lyapunov_Spec.txt')
current_dir = pwd;
if ~exist(sprintf('%s/Rossler_Info', current_dir), 'dir')
    mkdir(sprintf('%s/Rossler_Info', current_dir))
end

info = {'discarded', time_window+sum2; 'MLE', LYA(end,1); 'niter', k;'std_threshold', std_threshold; ...
  'coupling strength', coupling_strength; 'time window', time_window; 'max_discard', max_discard};
writecell(info, sprintf('Rossler_Info/info_%.3f.txt', coupling_strength))
writematrix(tvalues, sprintf('Rossler_Info/tvals_%.3f.txt', coupling_strength))
writematrix(LYA(:,1), sprintf('Rossler_Info/MLE_%.3f.txt', coupling_strength))

discard = time_window+sum2; %length(LYA(:,1)) - k;
figure;
loglog(tvalues, LYA(:,1))
xlim([0 length(LYA(:,1))])
ylim([-abs(min(LYA(:,1))), 1.5*max(LYA(:,1))])
xline(discard*stept, 'r--')
yline(0, 'r-')
hold on
patch([0 log(discard) log(discard) 0], [-abs(log(min(LYA(:,1)))) -abs(log(min(LYA(:,1)))) 1.5*log(max(LYA(:,1))) 1.5*log(max(LYA(:,1)))],'red', 'FaceAlpha', 0.3)
t = text((log(discard))/2,-abs(log(min(LYA(:,1))))/2, 'Transient Period');
t.Rotation = 90;
hold off
saveas(gca, sprintf('Rossler_Info/MLE_over_Time_%.3f.png', coupling_strength))
close all
end


%% %%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%

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

function output= rossler_ext(t, X, Adj,w, a, f, c, K)

M = size(Adj,2);
x= X(1:M); y = X(M+1:2*M); z=X(2*M+1:3*M);
output = zeros(3*M+9*M^2,1);

 for i=1:M
    sum1 = 0;
        for j=1:M
            sum1 = sum1 + Adj(i,j)*(x(j)-x(i));
        end
    output(i) =  -w(i)*y(i)-z(i)+ K*sum1;
    output(M+i) = w(i)*x(i) + a*y(i);
    output(2*M+i) = f + z(i) * (x(i)-c);
 end

 %%%% Partial derivatives
x_dx = zeros(M,M); x_dy = zeros(M,M); x_dz = zeros(M,M);
y_dx = zeros(M,M); y_dy = zeros(M,M); y_dz = zeros(M,M);
z_dx = zeros(M,M); z_dy = zeros(M,M); z_dz = zeros(M,M);

for i=1:M
    for j=1:M
        if i ==j
        x_dx(i,j) = -K*sum(Adj(i,:));
        else
            x_dx(i,j) = K*Adj(i,j);
        end
    end
end
 for i=1:M
     x_dy(i,i)=-w(i);
 end
for i=1:M
     x_dz(i,i)=-1;
end
for i=1:M
     y_dx(i,i) = w(i);
end
for i=1:M
     y_dy(i,i)=a;
end
for i=1:M
     z_dx(i,i)= z(i);
end
for i=1:M
     z_dz(i,i)= x(i)-c;
end

Jacob_Lorenz = [x_dx, x_dy, x_dz; y_dx, y_dy, y_dz; z_dx, z_dy, z_dz];

Y = zeros(3*M,3*M);
for i=1:3*M
Y(:,i) = X(i*3*M+1:(i+1)*3*M);
end

output(3*M+1:3*M+(3*M)^2)= Jacob_Lorenz*Y;

end











