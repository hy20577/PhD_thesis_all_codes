

function [tvalues, LYA, Orbits] = FHN_generator(eps, a, sigma, r, phi, x0, t0, tf, k, time_window, std_threshold, max_discard)

%%% ode45 function is used to solve ode equation, needed to revise for
%%% other solution.

% ext_odefun - extended ode system needed as input
% x0: initial condition for the ode system
% t0: initial time
% tf: final time
% k: number of iteration
% How many data you want to discard in computation of LE as 

if max_discard < time_window
    max_discard = time_window;
    warning('max_discard = time_window set up automatically')
end


n = length(x0); % dimension of extended ode function.
M = n/2; 
% M = M(M>0); % Number of oscillators
R = r*M;
Adj = ringAdjMat(M,R); % generating the regular network

B = [cos(phi) sin(phi); 
    -sin(phi) cos(phi)];

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
        [T, Y] = ode45(@(t,y) FHN_ext(t,y, eps, a, sigma, B, Adj, M, R), tspan, y0, options);

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
        [T, Y] = ode45(@(t,y) FHN_ext(t,y, eps, a, sigma, B, Adj, M, R), tspan, y0, options);
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
        [T, Y] = ode45(@(t,y) FHN_ext(t,y, eps, a, sigma, B, Adj, M, R), tspan, y0, options);
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
% 
% currentpath = pwd;
% folder = sprintf('%s/Lorenz_Info', currentpath);
% if ~exist(folder, 'dir')
%     mkdir(folder)
% end
% 
% info = {'discarded', time_window+sum2; 'MLE', LYA(end,1); 'niter', k;'std_threshold', std_threshold; ...
%   'coupling strength', coupling_strength; 'time window', time_window; 'max_discard', max_discard};
% writecell(info, sprintf('Lorenz_Info/info_%.2f.txt', coupling_strength))
% writematrix(tvalues, sprintf('Lorenz_Info/tvals_%.2f.txt', coupling_strength))
% writematrix(LYA(:,1), sprintf('Lorenz_Info/MLE_%.2f.txt', coupling_strength))
% 
% discard = time_window+sum2; %length(LYA(:,1)) - k;
% figure;
% loglog(tvalues, LYA(:,1))
% xlim([0 length(LYA(:,1))])
% ylim([-abs(min(LYA(:,1))), 1.5*max(LYA(:,1))])
% xline(discard*stept, 'r--')
% yline(0, 'r-')
% hold on
% patch([0 discard discard 0], [-abs(min(LYA(:,1))) -abs(min(LYA(:,1))) 1.5*max(LYA(:,1)) 1.5*max(LYA(:,1))],'red', 'FaceAlpha', 0.3)
% t = text((discard)/2,-abs(min(LYA(:,1)))/2, 'Transient Period');
% t.Rotation = 90;
% hold off
% saveas(gca, sprintf('Lorenz_Info/MLE_over_Time_%.2f.png', coupling_strength))
% close all
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


function f=FHN_ext(t,X,eps, a, sigma, B, Adj, M, R)
% sigma coupling strength

u=X(1:M); v=X(M+1:2*M);
f=zeros(2*M+4*M^2,1);

if length(a) == 1
    a = ones(M,1).*a;
end

        for k=1:M
            sum1 = 0;
            sum2 = 0;
                for j=1:M
                    sum1 = sum1 + Adj(k,j)*(B(1,1)*(u(j)-u(k))+B(1,2)*v(j)-v(k));
                    sum2 = sum2 + Adj(k,j)*(B(2,1)*(u(j)-u(k))+B(2,2)*v(j)-v(k));
                end
            f(k) = u(k) - (u(k)^3)/3 - v(k)+ (sigma/(2*R))*sum1;
            f(M+k) = u(k)+a(k) + (sigma/(2*R))*sum2;
        end

u_du = zeros(M,M);
    for i=1:M
        for j=1:M
            if i == j
                u_du(i,j) = (1/eps)*(1-u(i)^2-(sigma/2*R)*B(1,1)*sum(Adj(i,:)));
            else
                u_du(i,j) = (1/eps)*(-B(1,1));
            end
        end
    end

 u_dv = zeros(M,M);
     for i=1:M
        for j=1:M
            if i == j
                u_dv(i,j) = (1/eps)*(-1-(sigma/2*R)*B(1,2)*sum(Adj(i,:)));
            else
                u_dv(i,j) = (1/eps)*(-B(1,2));
            end
        end
    end
 v_du = zeros(M,M);
    for i=1:M
        for j=1:M
            if i == j
                v_du(i,j) = 1-(sigma/2*R)*B(2,1)*sum(Adj(i,:));
            else
                v_du(i,j) = -B(2,1);
            end
        end
    end
 v_dv = zeros(M,M);
     for i=1:M
        for j=1:M
            if i == j
                v_dv(i,j) = -(sigma/2*R)*B(2,2)*sum(Adj(i,:));
            else
                v_dv(i,j) = -B(2,2);
            end
        end
     end

    Jacob_Lorenz = [u_du, u_dv; v_du v_dv];

Y = zeros(2*M,2*M);
for i=1:2*M
Y(:,i) = X(i*2*M+1:(i+1)*2*M);
end

f(2*M+1:2*M+(2*M)^2)= Jacob_Lorenz*Y;

end




