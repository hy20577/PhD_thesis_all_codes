function [tvalues, LYA, Orbits] = HR_generator(x0, t0, tf, k, coupling_strength, Adj, time_window, std_threshold, max_discard)

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
            [T, Y] = ode45(@(t,y) HR_ext(t,y, Adj, coupling_strength), tspan, y0, options);
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
            [T, Y] = ode45(@(t,y) HR_ext(t,y, Adj, coupling_strength), tspan, y0, options);
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
    
         if sum2 > max_discard-time_window 
            warning(sprintf('std_threshold is too small to stabilize; first %d discarded', max_discard))
            break
         end
     end

     for ITERLAP=1:k   %Main loop to compute Lyapunov
        
        tspan= [t t+stept];
        options = odeset('RelTol',1e-8);
        [T, Y] = ode45(@(t,y) HR_ext(t,y, Adj, coupling_strength), tspan, y0, options);
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


currentpath = pwd;
folder = sprintf('%s/HR_Info', currentpath);
if ~exist(folder, 'dir')
    mkdir(folder)
end

info = {'discarded', time_window+sum2; 'MLE', LYA(end,1); 'niter', k;'std_threshold', std_threshold; ...
  'coupling strength', coupling_strength; 'time window', time_window; 'max_discard', max_discard};
writecell(info, sprintf('HR_Info/info_%.2f.txt', coupling_strength))
writematrix(tvalues, sprintf('HR_Info/tvals_%.2f.txt', coupling_strength))
writematrix(LYA(:,1), sprintf('HR_Info/MLE_%.2f.txt', coupling_strength))

discard = time_window+sum2; %length(LYA(:,1)) - k;
figure;
loglog(tvalues, LYA(:,1))
xlim([0 length(LYA(:,1))])
ylim([-abs(min(LYA(:,1))), 1.5*max(LYA(:,1))])
xline(discard*stept, 'r--')
yline(0, 'r-')
hold on
patch([0 discard discard 0], [-abs(min(LYA(:,1))) -abs(min(LYA(:,1))) 1.5*max(LYA(:,1)) 1.5*max(LYA(:,1))],'red', 'FaceAlpha', 0.3)
t = text((discard)/2,-abs(min(LYA(:,1)))/2, 'Transient Period');
t.Rotation = 90;
hold off
saveas(gca, sprintf('HR_Info/MLE_over_Time_%.2f.png', coupling_strength))
close all
end







%%%% FUNCTIONS   %%%%%%%%%

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
function f=HR_ext(t,y, A, g_l)
% We need the extended function to evolve deviation vectors in odesolver
% simultaneously. Detail in Chaos chapter 9. 

%%% Get ready inputs
M=size(A,1);
D = zeros(M,M);   %%% diagonal matrix of degree nodes in eq.9

    for j=1:M
    D(j,j) = sum(A(j,:));    
    end

C = D - A;   %%% Laplacian Matrix
a=1;
b=3;
d=5;
s=4;
h = 0.005;

f= zeros(12,1);

%%% Ode part
 p_dot = zeros(1,M);
 q_dot = zeros(1,M);
 n_dot = zeros(1,M);
    for i=1:M

        sum1 = 0;
            for l=1:M
            sum1 = sum1 + C(i,l)*y(l);
            end
        
        p_dot(i)  =  y(i+M) - y(i)^3 + 3*y(i)^2 - y(i+2*M) + 3.25 - g_l*sum1; 
        q_dot(i)  =  1 - 5*y(i)^2 - y(i+M); 
        n_dot(i)  =  0.02*(y(i)+1.6) - 0.005*y(i+2*M);
    end
    f(1:3*M) = [p_dot q_dot  n_dot];
    
% Jacobian Matrix

p_dp = zeros(M, M);
for j=1:M
    for i=1:M
        if i == j
         p_dp(i,j) = -3*a*y(i)^2 + b*y(i) - g_l* C(i,i);
        else
         p_dp(i,j) = - g_l*C(i,j);
        end
    end
end

q_dp = zeros(M,M);
    for i=1:M
        q_dp(i,i) = -2*d*y(i);
    end

n_dp = zeros(M,M);
    for i= 1:M
        n_dp(i,i)= h*s;
    end
p_dq =eye(M);
q_dq = -eye(M);
n_dq = zeros(M,M);
p_dn = -eye(M);
q_dn = zeros(M,M);
n_dn = -h*eye(M);

Jacob_HR = [p_dp p_dq p_dn; q_dp q_dq q_dn; n_dp n_dq n_dn]; % 36x36 Jacobian
Y = zeros(3*M,3*M);

for i=1:3*M
Y(:,i) = y(i*3*M+1:(i+1)*3*M);
end

f(3*M+1:3*M+(3*M)^2)= Jacob_HR*Y;
end