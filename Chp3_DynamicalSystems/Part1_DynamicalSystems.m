
clc; clear; close;

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp1 = strjoin(splittedstr(1:end-1), '/');
addpath(genpath(sprintf('%s/Functions', pp1)))


fig_folder = sprintf('%s/Figures',current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end

%% 1. Undamped single pendulum phase portrait


x_array = linspace(-2*pi,2*pi,15);
y_array = linspace(-10,10,15);
starting_pts = zeros(length(x_array)*length(y_array),2);
counter =0;
figure
for i=1:length(x_array)
    for j=1:length(y_array)
        counter= counter +1;
        [t,y] = ode45(@sing_pend,[0 10],[x_array(i); y_array(j)]);
        hold on
        plot(y(:,1),y(:,2),'k-')
        starting_pts(counter,1) = x_array(i);
        starting_pts(counter,2) = y_array(j);
    end
end

xlim([-2*pi-.1, 2*pi+.1])
ylim([-11,10])
yline(0)

hold on
quiver(starting_pts(:,1), starting_pts(:,2), starting_pts(:,2), -9.87.*sin(starting_pts(:,1)),'r', 'LineWidth', 2, 'AutoScale','on')

xt=[-2*pi, -pi, 0 pi 2*pi].';
yt=[-10 -5 0 5 10].';
% set(gca,'XTick',xt,'XTickLabel','', ...
%           'YTick',yt,'YTickLabel',num2str(yt,'%d'))
% text(xt,-11*ones(size(xt)),{'-2\pi', '-\pi', '0', '\pi', '2\pi'}, ...
%           'horizontal','center','vertical','top')

ax = gca;
ax.XTick = xt;
ax.XTickLabel = {'-2\pi', '-\pi', '0', '\pi', '2\pi'};
ax.YTick = yt;
ax.YTickLabel = arrayfun(@num2str, yt, 'UniformOutput', false);
ax.TickDir = "out";
xlabel('$\theta$', 'Interpreter','latex', 'FontSize',30)
ylabel('$\frac{d\theta}{dt}$', 'interpreter', 'latex', 'FontSize',30)
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/single_pendulum_phasePortrait.fig', fig_folder))
saveas(gcf, sprintf('%s/single_pendulum_phasePortrait', fig_folder), 'epsc')

%% 2. Population function of time

t=0:0.001:10;

x=2.^t;
y=2*t;

figure;
subplot(1,2,1)
plot(t,x, 'LineWidth',2);
xlabel('$t$','Interpreter','latex', 'FontSize', 30);
ylabel('$f(t)$','Interpreter','latex','FontSize', 30);
pbaspect([1 1 1])
title('(a)', 'Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

subplot(1,2,2)
plot(t,y, 'LineWidth',2);
xlabel('$X_{n}$','Interpreter','latex','FontSize', 30)
ylabel('$X_{n+1}$','Interpreter','latex','FontSize', 30)
pbaspect([1 1 1])
title('(b)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/PopulationFunction.fig', fig_folder))
saveas(gcf, sprintf('%s/PopulationFunction', fig_folder), 'epsc')

%% 3. Cobweb Plot   %%%%%%%%%%%%

r= 2;
% k=2;
a=0;
b=1;
x0= 0.01;
n=50;

x = x0;
orbits = zeros(n,1);
for iter = 1:n
orbits(iter) =x;
x = logis(x,r);
end

figure
subplot(2,2,1)
plot(orbits, 'Linewidth',2)
ylim([0,1])
title('(a)', 'FontSize',30, 'Interpreter','latex')
xlabel('$n$', 'FontSize',30, 'Interpreter','latex')
ylabel('$x_{n}$', 'FontSize',30, 'Interpreter','latex')
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;


subplot(2,2,2)
cobweb(@logis,a,b,x0,n, r)
title('(b)' ,'Interpreter','latex','FontSize',30)
xlabel('$x_{n}$','Interpreter','latex','FontSize',30)
ylabel('$x_{n+1}$','Interpreter','latex','FontSize',30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

%  2-Periodic Orbits

k=1;
r=3.4; 
% n=100; 
x0=0.4; % initial x
orbits= zeros(n,1);

x = x0;
for i=1:n
orbits(i) = logis(x,r); 
x = orbits(i);
end

subplot(2,2,3)
plot(1:n,orbits,'*-', 'MarkerSize',10)
title('(c)','Interpreter','latex','Fontsize',30) %2-Periodic Orbit of Logistic Map with $r$ =3.4
xlabel('$n$', 'Interpreter','latex','Fontsize',30)
ylabel('$x_{n}$', 'FontSize',30, 'Interpreter','latex')
ylim([0,1])
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

subplot(2,2,4)
cobweb(@logis, 0, 1, x0, n, r)
title('(d)','Interpreter','latex','FontSize',30)
xlabel('$x_{n}$','Interpreter','latex','FontSize',30)
ylabel('$x_{n+1}$','Interpreter','latex','FontSize',30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

%%
saveas(gcf, sprintf('%s/Cobweb_periodic', fig_folder))
saveas(gcf, sprintf('%s/Cobweb_periodic', fig_folder), 'epsc')

%% 4. Bifurcation of Logistic Map

n=1e4;
a=linspace(1,4,n);
t=1e3;
B = zeros(t-900,n);


for j=1:n
x=rand(1,1);
A = zeros(t,1);
    for i=1:t
     A(i) = x;
     x = logis(x, a(j));
    end
    B(:,j) = A(901:end);
end
%
figure
for k=1:n
plot(repmat(a(k),t-900,1),B(:,k),'k.', 'MarkerSize', 1)
hold on
end

ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;
xlabel('r', 'FontSize',30, 'Interpreter','latex')

%%
saveas(gcf, sprintf('%s/LogisticBifurcation.fig', fig_folder))
saveas(gcf, sprintf('%s/LogisticBifurcation', fig_folder), 'epsc')

%% 5. Linear & Nonlinear Transformation

M=0;
r=1;
n=1e3;
t=linspace(0,2*pi,n);
x=M+r*cos(t);
y=M+r*sin(t);


%%%% Linear Transformation
T = [1 ,2; 2,1];
B = T*[x; y];
C = T *B;
%%% Nonlinear Transformation
%
for i=1:n
K(:,i) = [1, 2*x(i);2, 1] * [x(i);y(i)];
end

plot(K(1,:),K(2,:))


tiledlayout(1,3)
nexttile   %% Linear Transformation
plot(x,y)
hold on
plot(B(1,:), B(2,:))
pbaspect([1 1 1])
hold on 
plot(C(1,:), C(2,:))
title('(a)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;
xlabel('$x$', 'FontSize',30, 'Interpreter','latex')
ylabel('$y$', 'FontSize',30, 'Interpreter','latex')


nexttile
plot(x,y)
hold on
plot(K(1,:),K(2,:))
pbaspect([1 1 1])
title('(b)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;
xlabel('$x$', 'FontSize',30, 'Interpreter','latex')
ylabel('$y$', 'FontSize',30, 'Interpreter','latex')

M=0;
r=0.1;
n=1e3;
t=linspace(0,2*pi,n);
x1=M+r*cos(t);
y1=M+r*sin(t);


for i=1:n
L(:,i) = [1, 2*x1(i); 2, 1] * [x1(i);y1(i)];
end

nexttile
plot(x1,y1)
xlim([-1 1])
ylim([-1 1])
hold on
plot(L(1,:),L(2,:))
pbaspect([1 1 1])
title('(c)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;
xlabel('$x$', 'FontSize',30, 'Interpreter','latex')
ylabel('$y$', 'FontSize',30, 'Interpreter','latex')
%%
saveas(gcf, sprintf('%s/Transformations.fig', fig_folder))
saveas(gcf, sprintf('%s/Transformations', fig_folder), 'epsc')

%% 6. 2D Lyapunov- evolution of variational vectors 

M=0;
r=1;
n=1e3;
t=linspace(0,2*pi,n);
x=M+r*cos(t);
y=M+r*sin(t);
for i=1:n
J = [4, -8/5; 4, -16/5];
trans(i,:) = J*[x(i); y(i)];
end



p0 = [0 0];
p1 = [1 0];
p2 = [0 1];


% Jacobian Matrix at initial point x=0.2, y=0.4.

J = [4, -8/5; 4, -16/5];

p1_iter = J*p1';
p2_iter = J*p2';


A=gramschmidt([p1_iter, p2_iter]);
p1_ortho = A(:,1);
p2_ortho = A(:,2);

p1_norm = p1_ortho/norm(p1_ortho);
p2_norm = p2_ortho/norm(p2_ortho);

figure(1);

subplot(2,2,1)
vectarrow(p0,p1)
hold on
vectarrow(p0,p2)
hold on
plot(x,y)
xlim([-6 6])
ylim([-6 6])
grid on
axis equal
title('(a)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

subplot(2,2,2)
vectarrow(p0,p1_iter);
hold on 
vectarrow(p0,p2_iter)
hold on 
plot(trans(:,1),trans(:,2))
xlim([-6 6])
ylim([-6 6])
axis equal
title('(b)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

subplot(2,2,3)
vectarrow(p0,p1_ortho);
hold on 
vectarrow(p0,p2_ortho)
hold on 
plot(trans(:,1),trans(:,2))
xlim([-6 6])
ylim([-6 6])
axis equal
title('(c)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;

subplot(2,2,4)
vectarrow(p0,p1_norm);
hold on 
vectarrow(p0,p2_norm)
hold on
plot(x,y)
xlim([-6 6])
ylim([-6 6])
axis equal
title('(d)','Interpreter','latex', 'FontSize', 30)
ax = gca;
ax.TickDir = "out";
ax.FontSize = 30;
%%
saveas(gcf, sprintf('%s/VariationVectors.fig', fig_folder))
saveas(gcf, sprintf('%s/VariationVectors', fig_folder), 'epsc')



%%   %%%%%%%%%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%


function cobweb(f, a, b, x0, n, varargin)
% f ---> function
% (a,b) ---> interval to show graph
% x0 ---> initial x value
% n ---> Number of iteration
x = linspace(a,b,1e5);
y = f(x, varargin{1});
plot(x, y, 'k', 'Linewidth',2);
hold on
plot(x,x, 'r', 'Linewidth',2);
line([x0,x0],[0,f(x0, varargin{1})], 'Linewidth',2);

    for i=1:n
        line([x0,f(x0, varargin{1})], [f(x0, varargin{1}),f(x0, varargin{1})]);
        x0 = f(x0, varargin{1});
        line([x0,x0],[x0,f(x0, varargin{1})], 'Linewidth',2)
    end
xlabel('$x_{n}$','Interpreter','latex','FontSize',20)
ylabel('$x_{n+1}$','Interpreter','latex','FontSize',20)
pbaspect([1 1 1])
end

function l=logis(x,r)
% global r
l=r.*x.*(1-x);
end

% function comp_l = comp2(x,r)
% %%% nth order Composite Logistic Function.
% n=2;
%     for i=1:n
%         x_val = r.*x.*(1-x);
%         x =x_val;
%     end
% comp_l = x_val;
% end
 

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

function vectarrow(p0,p1)
%Arrowline 3-D vector plot.
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1. The function can plot both 2D and 3D vector with arrow
%   depending on the dimension of the input
%
%   Example:
%       3D vector
%       p0 = [1 2 3];   % Coordinate of the first point p0
%       p1 = [4 5 6];   % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%   See also Vectline
%   Rentian Xiong 4-18-05
%   $Revision: 1.0
  if max(size(p0))==3
      if max(size(p1))==3
          x0 = p0(1);
          y0 = p0(2);
          z0 = p0(3);
          x1 = p1(1);
          y1 = p1(2);
          z1 = p1(3);
          plot3([x0;x1],[y0;y1],[z0;z1]);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.1;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
          
          hold on
          plot3(hu(:),hv(:),hw(:))  % Plot arrow head
          grid on
          % xlabel('x')
          % ylabel('y')
          % zlabel('z')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  elseif max(size(p0))==2
      if max(size(p1))==2
          x0 = p0(1);
          y0 = p0(2);
          x1 = p1(1);
          y1 = p1(2);
          plot([x0;x1],[y0;y1]);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.1;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          
          hold on
          plot(hu(:),hv(:))  % Plot arrow head
          grid on
          % xlabel('x')
          % ylabel('y')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  else
      error('this function only accepts 2D or 3D vector')
  end
end

