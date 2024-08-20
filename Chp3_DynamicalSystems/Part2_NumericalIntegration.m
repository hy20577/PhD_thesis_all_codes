
%% Numerical Integration Methods

clc; clear; close;

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp1 = strjoin(splittedstr(1:end-1), '/');
addpath(genpath(sprintf('%s/Functions', pp1)))


fig_folder = sprintf('%s/Figures',current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end


%%
%%%%%%%%%%%%%%%%%%%%%%   PART 1: ESSENTIALS for ALL MODELS  %%%%%%%%%%%%%%%%%%%%
%% Parameters

%%%%%   Constant parameters of HR Neuron Model

a=1;
b=3;
c=1;
d=5;
s=4;
p_0 = -1.6;
I_ext = 3.25;
h = 0.005;

%%%%% Parameters to set up

g_l= 0 ; %% Coupling Strength
t0=0;   %% Initial time point
tf=1000; %% Final time point

%%%%%%%%%%%%%% Laplacian Matrix %%%%%%%%%

A = readmatrix(sprintf('%s/Functions/12nodes.txt', pp1));
M=size(A,1);
K = zeros(M,M);   %%% diagonal matrix of degree nodes in eq.9

for j=1:M
K(j,j) = sum(A(j,:));    
end

C = K - A;   %%% Laplacian Matrix

% Initial Conditions

etap = 0.5*rand(M,1);
etaq = 0.5*rand(M,1);
etan = 0.5*rand(M,1);

p_init = zeros(1,M);
q_init = zeros(1,M);
n_init = zeros(1,M);


for i=1:M
        p_init(i)   = -1.30784489 + etap(i);
        q_init(i)   = -7.32183132 + etaq(i);      %% Setting initial conditions
        n_init(i)   =  3.35299859 + etan(i);
end


%%   PART 2: MODEL APPLICATIONS  %%%%%%%%%%%%%%%%%%%%

%% 2.1.EULER Approximation

%%%%%%% Model Parameters %%%%%%%%%%%%
starting_time = tic();
k = 10^7; %% number of iterations
step_size = (tf-t0)/k;

%%%%%%%% Initial Conditions %%%%%%%%

E_p = zeros(k,M);
E_q = zeros(k,M);
E_n = zeros(k,M);
tspan_E = zeros(k,1); % this is the time vector

E_p(1,:) = p_init;
E_q(1,:) = q_init;
E_n(1,:) = n_init;

%%%%%%%% Model Application  %%%%%%%

for j=1:k-1   % time loop
        E_p(j+1,:) = E_p(j,:) + step_size * f(E_p(j,:),E_q(j,:),E_n(j,:),C,g_l,a,b,I_ext,M);
        E_q(j+1,:) = E_q(j,:) + step_size * g(E_p(j,:),E_q(j,:),c,d);
        E_n(j+1,:) = E_n(j,:) + step_size * z(E_p(j,:),E_n(j,:),h,s,p_0);
        tspan_E(j+1) = step_size*j;
end

sprintf('The running time of Euler approximation: %.5f secs', toc(starting_time))


%% 2.2. RUNGE-KUTTA APPROXIMATION

% 4th Order-RK model is used.

% Model Parameters
starting_time = tic();
k = 10^6; %% number of iterations
step_size = (tf-t0)/k;

%%%%%%%%%%% Initial Conditions %%%%%%%%%

p_RK = zeros(k,M);
q_RK = zeros(k,M);
n_RK = zeros(k,M);
tspan_RK = zeros(k,1);

p_RK(1,:) = p_init;
q_RK(1,:) = q_init;   %% setting same initial conditions
n_RK(1,:) = n_init;

%%%%%%%%%% RK4 application  %%%%%%%%%%

for j=1:k-1   % time loop
          p_1  = f(p_RK(j,:),q_RK(j,:),n_RK(j,:),C,g_l,a,b,I_ext,M);
          q_1  = g(p_RK(j,:),q_RK(j,:),c,d);
          n_1  = z(p_RK(j,:),n_RK(j,:),h,s,p_0);
         
    
          p_2  = f(p_RK(j,:)+(step_size/2)*p_1, q_RK(j,:)+(step_size/2)*q_1, n_RK(j,:)+(step_size/2)*n_1,C,g_l,a,b,I_ext,M);
          q_2  = g(p_RK(j,:)+(step_size/2)*p_1, q_RK(j,:)+(step_size/2)*q_1,c,d);
          n_2  = z(p_RK(j,:)+(step_size/2)*p_1, n_RK(j,:)+(step_size/2)*n_1,h,s,p_0);
        
               
          p_3  = f(p_RK(j,:)+(step_size/2)*p_2, q_RK(j,:)+(step_size/2)*q_2, n_RK(j,:)+(step_size/2)*n_2,C,g_l,a,b,I_ext,M);
          q_3  = g(p_RK(j,:)+(step_size/2)*p_2, q_RK(j,:)+(step_size/2)*q_2,c,d);
          n_3  = z(p_RK(j,:)+(step_size/2)*p_2, n_RK(j,:)+(step_size/2)*n_2,h,s,p_0);
        
                
          p_4  = f(p_RK(j,:)+(step_size)*p_3, q_RK(j,:)+(step_size)*q_3, n_RK(j,:)+(step_size)*n_3,C,g_l,a,b,I_ext,M);
          q_4  = g(p_RK(j,:)+(step_size)*p_3, q_RK(j,:)+(step_size)*q_3,c,d);
          n_4  = z(p_RK(j,:)+(step_size)*p_3, n_RK(j,:)+(step_size)*n_3,h,s,p_0);
        
          p_RK(j+1,:) = p_RK(j,:) + (step_size/6) * (p_1+2*p_2+2*p_3+p_4);
          q_RK(j+1,:) = q_RK(j,:) + (step_size/6) * (q_1+2*q_2+2*q_3+q_4);
          n_RK(j+1,:) = n_RK(j,:) + (step_size/6) * (n_1+2*n_2+2*n_3+n_4);
          
          tspan_RK(j+1) = step_size*j;
end

sprintf('The running time of RK4 approximation: %.5f secs', toc(starting_time))

%% 2.3.Multistep: Adams-Bashford
% 4th order Adam Bashford based on previous 4 RK approximations.
starting_time = tic();

%%%%%%%%%% Initial Conditions  %%%%%%%%%
 
p_AB = zeros(k,M);
q_AB = zeros(k,M);
n_AB = zeros(k,M);

p_AB(1:4,:)=p_RK(1:4,:);
q_AB(1:4,:)=q_RK(1:4,:);   %%% First 4 time steps of AB is exactly same with
n_AB(1:4,:)=n_RK(1:4,:);    %% the first 4 time steps of RK!

tspan_AB = zeros(k,1);
tspan_AB(1:4,:) =  tspan_RK(1:4,:);

%%%%%%%%%% Model Application %%%%%%%%%%%%%%

for j=1:(k-4)
        p_n = f(p_AB(j,:), q_AB(j,:), n_RK(j,:),C,g_l,a,b,I_ext,M);
        p_n1 = f(p_AB(j+1,:), q_AB(j+1,:), n_AB(j+1,:),C,g_l,a,b,I_ext,M);
        p_n2 = f(p_AB(j+2,:), q_AB(j+2,:), n_AB(j+2,:),C,g_l,a,b,I_ext,M);
        p_n3 = f(p_AB(j+3,:), q_AB(j+3,:), n_AB(j+3,:),C,g_l,a,b,I_ext,M);
    
        q_n = g(p_AB(j,:), q_AB(j,:),c,d);
        q_n1 = g(p_AB(j+1,:), q_AB(j+1,:),c,d);
        q_n2 = g(p_AB(j+2,:), q_AB(j+2,:),c,d);
        q_n3 = g(p_AB(j+3,:), q_AB(j+3,:),c,d);
    
        n_n = z(p_AB(j,:) , n_AB(j,:) ,h,s,p_0);
        n_n1 = z(p_AB(j+1,:), n_AB(j+1,:) ,h,s,p_0);
        n_n2 = z(p_AB(j+2,:), n_AB(j+2,:) ,h,s,p_0);
        n_n3 = z(p_AB(j+3,:)  , n_AB(j+3,:)   ,h,s,p_0);

        p_AB(j+4,:) = p_AB(j+3,:) + (step_size/24)*(55*p_n3 - 59*p_n2+37*p_n1-9*p_n);
        q_AB(j+4,:) = q_AB(j+3,:) + (step_size/24)*(55*q_n3 - 59*q_n2+37*q_n1-9*q_n);
        n_AB(j+4,:) = n_AB(j+3,:) + (step_size/24)*(55*n_n3 - 59*n_n2+37*n_n1-9*n_n);
     
        tspan_AB(j+4) = step_size*(j+3);
end

sprintf('The running time of Adams-Bashford approximation: %.5f secs', toc(starting_time))

% 2.4. Ode23
starting_time = tic();

%%%%%%%%%% Initial Value  %%%%%
    
    k = 10^6;  %% Same iteration number with RK and AB  
    y_ode23 = zeros(k,3*M);
    tspan_ode23 = zeros(k,1);
    
    y_ode23(1,:) = [p_init q_init n_init];
    
%%%%%%%%% Model Application %%%%%%   
for i =1:(k-1)
        
        tspan = [(i-1)/1000 (i)/1000];
        options = odeset('RelTol',1e-12);
        [t,y] = ode23(@(t,y) odefun(t,y,M,C,g_l), tspan, y_ode23(i,:), options);
        
        y_ode23(i+1,:) = y(end,:);
        tspan_ode23(i+1) = (t(end));
end
    
    ode23_p = y_ode23(:,1:M);          %%% write p,q,n values into separate matrices.
    ode23_q = y_ode23(:,(M+1):2*M);
    ode23_n = y_ode23(:,(2*M+1):3*M);    

sprintf('The running time of ODE23 approximation: %.5f secs', toc(starting_time))

% 2.5. Ode45
starting_time = tic();

%%%%%%%% Initial Value %%%%%%%%
    k = 10^6;  %% Same iteration number with RK and AB  
    y_ode45 = zeros(k,3*M);
    tspan_ode45 = zeros(k,1);
    
    y_ode45(1,:) = [p_init q_init n_init];
    
%%%%%% Model Application
    

    for i =1:(k-1)
        
        tspan = [(i-1)/1000 (i)/1000];
        options = odeset('RelTol',1e-8);
        [t1,y1] = ode45(@(t1,y1) odefun(t1,y1,M,C,g_l), tspan, y_ode45(i,:), options);
        
        y_ode45(i+1,:) = y1(end,:);
        tspan_ode45(i+1) = t1(end); 
    end
    
  
    ode45_p = y_ode45(:,1:M);
    ode45_q = y_ode45(:,(M+1):2*M);
    ode45_n = y_ode45(:,(2*M+1):3*M);
    
sprintf('The running time of ODE45 approximation: %.5f secs', toc(starting_time))


%%                       PART 3: ABSOLUTE RELATIVE ERRORS  %%%%%%%%%%%%%%%%%%%%
    % To compare the methods, we use the sum of absolute relative error for
% each variable p,q, and n of M nodes, which is denoted by 'SRE'
    % To investigate the error distribution with log-log graph, we also
% compute the relative error in each time step for each M node, which is 
% denoted by 'RE'

%% 3.1. Runge-Kutta versus Adam-Bashford 
 
RE_RK_AB_p = zeros(length(tspan_RK),M);
RE_RK_AB_q = zeros(length(tspan_RK),M);
RE_RK_AB_n = zeros(length(tspan_RK),M);

SRE_RK_AB_p = zeros(M,1);
SRE_RK_AB_q = zeros(M,1);
SRE_RK_AB_n = zeros(M,1);


for i=1:M
RE_RK_AB_p(:,i) = are(p_AB(:,i), p_RK(:,i)); % Abs. rel. Error (are) between RK and AB
RE_RK_AB_q(:,i) = are(q_AB(:,i), q_RK(:,i));
RE_RK_AB_n(:,i) = are(n_AB(:,i), n_RK(:,i));

SRE_RK_AB_p(i) = sum(RE_RK_AB_p(:,i));
SRE_RK_AB_q(i) = sum(RE_RK_AB_q(:,i));
SRE_RK_AB_n(i) = sum(RE_RK_AB_n(:,i));

end


%% 3.2. Runge-Kutta versus Euler   

RE_RK_E_p = zeros(length(tspan_RK),M);
RE_RK_E_q = zeros(length(tspan_RK),M);
RE_RK_E_n = zeros(length(tspan_RK),M);

SRE_RK_E_p = zeros(M,1);
SRE_RK_E_q = zeros(M,1);
SRE_RK_E_n = zeros(M,1);

for i=1:M
RE_RK_E_p(:,i) = are(E_p(1:10:end,i), p_RK(:,i)); % Abs. rel. Error (are) between RK and Euler
RE_RK_E_q(:,i) = are(E_q(1:10:end,i), q_RK(:,i));
RE_RK_E_n(:,i) = are(E_n(1:10:end,i), n_RK(:,i));

SRE_RK_E_p(i) = sum(RE_RK_E_p(:,i)); % Abs. rel. Error (are) between RK and Euler
SRE_RK_E_q(i) = sum(RE_RK_E_q(:,i));
SRE_RK_E_n(i) = sum(RE_RK_E_n(:,i));

end 

%% 3.3. Adam-Bashford versus Euler  %%%%%%%%%%%

RE_AB_E_p = zeros(length(tspan_AB),M);
RE_AB_E_q = zeros(length(tspan_AB),M);
RE_AB_E_n = zeros(length(tspan_AB),M);

SRE_AB_E_p = zeros(M,1);
SRE_AB_E_q = zeros(M,1);
SRE_AB_E_n = zeros(M,1);

for i=1:M
RE_AB_E_p(:,i) = are(p_AB(:,i),E_p(1:10:end,i)); % Abs. rel. Error (are) between RK and Euler
RE_AB_E_q(:,i) = are(q_AB(:,i),E_q(1:10:end,i));
RE_AB_E_n(:,i) = are(n_AB(:,i),E_n(1:10:end,i));

SRE_AB_E_p(i) = sum(RE_AB_E_p(:,i));
SRE_AB_E_q(i) = sum(RE_AB_E_q(:,i));
SRE_AB_E_n(i) = sum(RE_AB_E_n(:,i));
end


%% 3.4. Ode23 versus Ode45


RE_ode23_ode45_p = zeros(length(tspan_ode45),M); 
RE_ode23_ode45_q = zeros(length(tspan_ode45),M);
RE_ode23_ode45_n = zeros(length(tspan_ode45),M); 

SRE_ode23_ode45_p = zeros(M,1);
SRE_ode23_ode45_q = zeros(M,1);
SRE_ode23_ode45_n = zeros(M,1);

for j=1:M
RE_ode23_ode45_p(:,j) = are(ode23_p(:,j),ode45_p(:,j)); 
RE_ode23_ode45_q(:,j) = are(ode23_q(:,j),ode45_q(:,j)); 
RE_ode23_ode45_n(:,j) = are(ode23_n(:,j),ode45_n(:,j)); 

SRE_ode23_ode45_p(j) = sum(RE_ode23_ode45_p(:,j));
SRE_ode23_ode45_q(j) = sum(RE_ode23_ode45_q(:,j));
SRE_ode23_ode45_n(j) = sum(RE_ode23_ode45_n(:,j));
end


    % We observe RK and Adam-Bashford give similar results but different 
% than Euler, it is because Euler cannot keep accuracy up for long
% integration time. Therefore, we compare RK with ode23 and ode45.


%% 3.5. RK versus Ode23 

RE_ode23_RK_p = zeros(length(tspan_RK),M); 
RE_ode23_RK_q = zeros(length(tspan_RK),M);
RE_ode23_RK_n = zeros(length(tspan_RK),M); 

SRE_ode23_RK_p = zeros(M,1);
SRE_ode23_RK_q = zeros(M,1);
SRE_ode23_RK_n = zeros(M,1);

for j=1:M
RE_ode23_RK_p(:,j) = are(ode23_p(:,j),p_RK(:,j)); 
RE_ode23_RK_q(:,j) = are(ode23_q(:,j),q_RK(:,j)); 
RE_ode23_RK_n(:,j) = are(ode23_n(:,j),n_RK(:,j)); 

SRE_ode23_RK_p(j) = sum(RE_ode23_RK_p(:,j));
SRE_ode23_RK_q(j) = sum(RE_ode23_RK_q(:,j));
SRE_ode23_RK_n(j) = sum(RE_ode23_RK_n(:,j));
end


%% 3.6. RK versus Ode45 

RE_ode45_RK_p = zeros(length(tspan_RK),M); 
RE_ode45_RK_q = zeros(length(tspan_RK),M);
RE_ode45_RK_n = zeros(length(tspan_RK),M); 

SRE_ode45_RK_p = zeros(M,1);
SRE_ode45_RK_q = zeros(M,1);
SRE_ode45_RK_n = zeros(M,1);

for j=1:M
RE_ode45_RK_p(:,j) = are(ode45_p(:,j),p_RK(:,j)); 
RE_ode45_RK_q(:,j) = are(ode45_q(:,j),q_RK(:,j)); 
RE_ode45_RK_n(:,j) = are(ode45_n(:,j),n_RK(:,j)); 

SRE_ode45_RK_p(j) = sum(RE_ode45_RK_p(:,j));
SRE_ode45_RK_q(j) = sum(RE_ode45_RK_q(:,j));
SRE_ode45_RK_n(j) = sum(RE_ode45_RK_n(:,j));
end

%% 3.7. AB versus Ode45

RE_ode45_AB_p = zeros(length(tspan_AB),M); 
RE_ode45_AB_q = zeros(length(tspan_AB),M);
RE_ode45_AB_n = zeros(length(tspan_AB),M); 

SRE_ode45_AB_p = zeros(M,1);
SRE_ode45_AB_q = zeros(M,1);
SRE_ode45_AB_n = zeros(M,1);

for j=1:M
RE_ode45_AB_p(:,j) = are(ode45_p(:,j),p_AB(:,j)); 
RE_ode45_AB_q(:,j) = are(ode45_q(:,j),q_AB(:,j)); 
RE_ode45_AB_n(:,j) = are(ode45_n(:,j),n_AB(:,j)); 

SRE_ode45_AB_p(j) = sum(RE_ode45_AB_p(:,j));
SRE_ode45_AB_q(j) = sum(RE_ode45_AB_q(:,j));
SRE_ode45_AB_n(j) = sum(RE_ode45_AB_n(:,j));
end



%%

Euler_p = [tspan_E, E_p];
Trajectories_RK_AB_RK23_RK45 = [tspan_RK,p_RK(:,1), p_AB(:,1), ode23_p(:,1), ode45_p(:,1)];
SRE


          %%    PART 4: PRESENTING THE RESULTS  %%
                      %% 4.1. GRAPHS





figure
subplot(4,2,1)
plot(tspan_RK,p_RK(:,1))
hold on
plot(tspan_E, E_p(:,1))
l1 = legend('RK4','Euler','Interpreter','latex', 'FontSize', 20);
title('(a)','Interpreter','latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 20; 

subplot(4,2,2)
plot((tspan_RK),cumsum(RE_RK_E_p(:,1)), 'LineWidth',2)
title('(e)','Interpreter','latex', 'FontSize', 20)
ylim([0,3e6]); 
ax = gca;
ax.FontSize = 20; 

subplot(4,2,3)
plot(tspan_RK,p_RK(:,1))
hold on
plot(tspan_AB,p_AB(:,1))
l2 = legend('RK4','AB4','Interpreter','latex', 'FontSize', 20);
title('(b)','Interpreter','latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 20; 

subplot(4,2,4)          %%% Error Distribution for p1
plot((tspan_AB),cumsum(RE_RK_AB_p(:,1)), 'LineWidth',2)
title('(f)','Interpreter','latex', 'FontSize', 20)
ylim([0,5e3])
ax = gca;
ax.FontSize = 20; 

subplot(4,2,5)
plot(tspan_RK,p_RK(:,1))
hold on
plot(tspan_ode23, ode23_p(:,1))
l3 =legend('RK4','RK23','Interpreter','latex', 'FontSize', 20, 'Location', 'bestoutside');
title('(c)','interpreter', 'latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 20;

subplot(4,2,6)
plot(tspan_RK,cumsum(RE_ode23_RK_p(:,1)), 'LineWidth',2)
title('(g)','Interpreter','latex', 'FontSize', 20)
ylim([0,1])
ax = gca;
ax.FontSize = 20;

subplot(4,2,7)
plot(tspan_RK,p_RK(:,1))
hold on
plot(tspan_ode45, ode45_p(:,1))
l4 = legend('RK4','RK45','Interpreter','latex', 'FontSize', 20);
title('(d)','Interpreter','latex', 'FontSize', 20)
ax = gca;
ax.FontSize = 20;

subplot(4,2,8)
plot(tspan_RK,cumsum(RE_ode45_RK_p(:,1)), 'LineWidth',2)
title('(h)','Interpreter','latex', 'FontSize', 20)
ylim([0,1])
ax = gca;
ax.FontSize = 20;
%%
saveas(gcf, sprintf('%s/IntegrationMethodsComparison.fig', fig_folder))
saveas(gcf, sprintf('%s/IntegrationMethodsComparison', fig_folder), 'epsc')
%%
% %% 4.1.1 RK versus AB
% 
% figure(1)           %%% Graph of RK and AB for p1
% tiledlayout(1,2)
% nexttile
% plot(tspan_RK,p_RK(:,1))
% hold on
% plot(tspan_AB,p_AB(:,1))
% legend('Runge-Kutta','Adams-Bashford','Interpreter','latex', 'FontSize', 30)
% title('Runge-Kutta versus Adams-Bashford','Interpreter','latex', 'FontSize', 30)
% 
% nexttile          %%% Error Distribution for p1
% loglog((tspan_AB),(RE_RK_AB_p(:,1)))
% title('Error Distribution of RK and AB: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% %% 4.1.2 RK versus Euler
% figure(2)
% tiledlayout(1,2)
% nexttile
% plot(tspan_RK,p_RK(:,1))
% hold on
% plot(tspan_E, E_p(:,1))
% legend('RK4','Euler','Interpreter','latex', 'FontSize', 30)
% title('RK4 versus Euler','Interpreter','latex', 'FontSize', 30)
% 
% nexttile
% loglog((tspan_RK),(RE_RK_E_p(:,1)))
% title('Error Distribution of RK4 and Euler: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% 
% %% 4.1.3 AB versus Euler
% 
% figure(3)
% tiledlayout(1,2)
% nexttile
% plot(tspan_AB,p_AB(:,1))
% hold on
% plot(tspan_E, E_p(:,1))
% legend('AB','Euler','Interpreter','latex', 'FontSize', 30)
% title('AB versus Euler','Interpreter','latex', 'FontSize', 30)
% 
% nexttile
% loglog(tspan_AB,RE_AB_E_p(:,1))
% title('Error Distribution of AB & Euler: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% 
% 
% %% 4.1.4 Ode23 versus Ode45
% figure(4)
% tiledlayout(1,2)
% nexttile
% plot(tspan_ode23,ode23_p(:,1))
% hold on
% plot(tspan_ode45, ode45_p(:,1))
% legend('ODE23','ODE45','Interpreter','latex', 'FontSize', 30)
% title('ODE23 versus ODE45','Interpreter','latex', 'FontSize', 30)
% 
% nexttile
% loglog(tspan_ode23, RE_ode23_ode45_p(:,1))
% title('Error Distribution of ODE23 and ODE45: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% 
% 
% %% 4.1.5 RK versus Ode23
% 
% figure(5)
% tiledlayout(1,2)
% nexttile
% plot(tspan_RK,p_RK(:,1))
% hold on
% plot(tspan_ode23, ode23_p(:,1))
% legend('RK','ODE23','Interpreter','latex', 'FontSize', 30)
% title('RK versus ODE23','Interpreter','latex', 'FontSize', 30)
% 
% nexttile
% loglog(tspan_RK,RE_ode23_RK_p(:,1))
% title('Error Distribution of RK and ODE23: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% 
% 
% %% 4.1.6 RK versus Ode45
% 
% figure(6)
% tiledlayout(1,2)
% nexttile
% plot(tspan_RK,p_RK(:,1))
% hold on
% plot(tspan_ode45, ode45_p(:,1))
% legend('RK','ODE45','Interpreter','latex', 'FontSize', 30)
% title('RK versus ODE45','Interpreter','latex', 'FontSize', 30)
% 
% nexttile
% loglog(tspan_RK,RE_ode45_RK_p(:,1))
% title('Error Distribution of RK and ODE45: Log-log graph','Interpreter','latex', 'FontSize', 30)
% 
% %% 4.1.7 AB versus Ode45
% figure(13)
% plot(tspan_AB,p_AB(:,1))
% hold on
% plot(tspan_ode45, ode45_p(:,1))
% legend('AB','ODE45')
% title('AB versus ODE45')
% 
% figure(14)
% loglog(tspan_AB,RE_ode45_AB_p(:,1))
% title('Error Distribution of AB & ODE45: Log-log graph')



                       %% 4.2 SUM of ERROR TABLES

%% 4.2.1 Runge-Kutta versus Adams-Bashford

SRE_RK_AB =[ SRE_RK_AB_p SRE_RK_AB_q SRE_RK_AB_n];

%% 4.2.2 Runge-Kutta versus Euler

SRE_RK_E =[ SRE_RK_E_p SRE_RK_E_q SRE_RK_E_n];

%% 4.2.3 Adams-Bashford versus Euler

SRE_AB_E =[ SRE_AB_E_p SRE_AB_E_q SRE_AB_E_n];

%% 4.2.4. Ode23 versus Ode45

SRE_ode23_ode45 =[ SRE_ode23_ode45_p SRE_ode23_ode45_q SRE_ode23_ode45_n];

%% 4.2.5 RK versus ode23

SRE_ode23_RK =[ SRE_ode23_RK_p SRE_ode23_RK_q SRE_ode23_RK_n];


%% 4.2.6 RK versus ode45

SRE_ode45_RK = [SRE_ode45_RK_p SRE_ode45_RK_q SRE_ode45_RK_n];


%% 4.2.7 AB versus Ode45

SRE_ode45_AB = [SRE_ode45_AB_p SRE_ode45_AB_q SRE_ode45_AB_n];






%% FUNCTIONS


%% HR model of neuronal bursting using 3-coupled first order differential equations

function pdot = f(p,q,n,C,g_l,a,b,I_ext,M)
    pdot= zeros(1,M);
    for i=1:M
        sum1 = 0;
        for l=1:M
            sum1 = sum1 + C(i,l)*p(l);
        end
        pdot(i)=q(i)-a*p(i)^3+b*p(i)^2-n(i)+ I_ext - g_l*sum1;
    end
end

function qdot=g(p,q,c,d)
qdot=c-d*p.^2-q;
end

function ndot=z(p,n,h,s,p_0)
ndot=h*(s*(p-p_0)-n);
end

%%%%%%%%%%%% Function for Ode23 & ode45

function dydt = odefun(~,y,M,C,g_l)
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
    dydt = [p_dot q_dot  n_dot]';
end



function abs_rel_err = are(a,b)
abs_rel_err = abs((a-b)./b);
end

%%%%%%%%%%%%%%%%%%%   Functions for SdeSolver  %%%%%%%%%%%%%%%%%%%  

function dr = drift(t,x, A, w)
% A is the adjaceny matrix
% theta is phases of oscilators
% K is the coupling strength
% global A w

K = 4;
N = size(A,2);
dr = zeros(N,1);

    for i=1:N
        sum1 =0; 
        for j=1:N
            sum1 = sum1+ A(i,j) * sin(x(j)-x(i));  
        end
        dr(i) = w(i) + (K/N)*sum1; 
    end
end
% 
function diffusion = diff(t,x, A, D)
% global A D
N = length(A(1,:));
diffusion= D*eye(N);
end

