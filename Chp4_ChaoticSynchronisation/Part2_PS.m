
%% The code produces the results presented in chapter 4: Synchronisation


close, clear, clc

current_path = pwd;
splittedstr = strsplit(current_path, '/');
pp = strjoin(splittedstr(1:end-1), '/'); 
addpath(genpath(sprintf(pp, '/Functions')))
% 

fig_folder = sprintf('%s/Figures', current_path); % create folder to save figures

if ~exist(fig_folder, 'dir')
    mkdir(sprintf('%s', fig_folder))
end

%% %% 1. Phase concept in periodic motion

% Monocomponent sin wave

f1 = 2; % frequency of the signal. f1 cycle in 1 sec. 
A1 = 2; % Amplitude of the signal 1.
Fs = 1e3; %sample rate, Fs, points in 1 sec.
t = 0:1/Fs:2-(1/Fs);

Dat1 = A1*sin(2*pi*f1*t);
%
figure
subplot(2,2,1:2)
plot(t,Dat1, 'LineWidth',2)
ylim([-4 4])
xlabel('$t$ (sec)', 'Interpreter', 'latex', 'Fontsize', 20)
ylabel('$x(t)$', 'Interpreter', 'latex', 'Fontsize', 20)

yline(0,'--k', 'LineWidth',2)
%
idx_pi_over2 = find(Dat1 == 2);
xpos = t(idx_pi_over2);
labels = {'\pi/2','5\pi/2', '9\pi/2', '13\pi/2'};
text(xpos-0.03, Dat1(idx_pi_over2)+0.5,labels, 'FontSize',30);
hold on
plot(xpos, Dat1(idx_pi_over2), 'or', 'LineWidth',2)
title('(a)', 'Interpreter', 'latex', 'Fontsize', 30)
ax = gca;
ax.FontSize = 30; 
xticks([0, 0.5, 1, 1.5, 2]);
yticks([-2, 0, 2]); 
%
subplot(2,2,3)
yyaxis left
plot(t,mod(2*pi*f1*t,2*pi), 'LineWidth',2)
xlabel('$t$ (sec)', 'Interpreter', 'latex', 'Fontsize', 20)
ylabel('$\theta(t)$', 'Interpreter', 'latex', 'Fontsize', 20)
yticks([0 pi 2*pi])
yticklabels({'0', '\pi', '2\pi'})
title('(b)', 'interpreter', 'latex', 'Fontsize', 30)
ylim([-7,7])
hold on

% plot(t,mod(2*pi*f1*t,2*pi) - mod(unwrap(hilbert_phases(Dat1)),2*pi) )

yyaxis right
plot(t , mod(unwrap(hilbert_phases(Dat1)),2*pi), 'LineWidth',2)
hold on 
plot(t , mod(unwrap(hilbert_phases(Dat1))+pi/2,2*pi), 'LineWidth',1,'LineStyle',':')
ylim([0,14])
yticks([0 pi 2*pi])
yticklabels({'0', '\pi', '2\pi'})
% ax.YTickLabel = {'0', '2\pi'};
ax = gca;
ax.FontSize = 30; 
ylabel('$\theta(t)$', 'Interpreter', 'latex', 'Fontsize', 20)
%
subplot(2,2,4)
plot(2*cos(2*pi*f1*t),2*sin(2*pi*f1*t), 'LineWidth',2)
xlim([-3, 3])
ylim([-3, 3])
xline(0, 'LineWidth',2); yline(0, 'LineWidth',2);
pbaspect([1 1 1])
grid on
hold on
plot([0,A1*cos(pi/4)], [0, A1*sin(pi/4)], 'LineWidth',2)
plot(A1*cos(pi/4),A1*sin(pi/4),'or', 'LineWidth',2)
title('(c)', 'interpreter', 'latex', 'Fontsize', 30)
ylabel('$A \sin \theta$', 'FontSize',20,'Interpreter', 'latex')
xlabel('$A \cos \theta$', 'FontSize',20,'Interpreter', 'latex')

plot(1/4*cos(1/16*pi*f1*t),1/4*sin(1/16*pi*f1*t), 'LineWidth',2)
text(0.5,0.3, '\pi/4', 'FontSize',20)
ax = gca;
ax.FontSize = 30; 

%%
saveas(gca, sprintf('%s/Phase_sinus_wave.fig', fig_folder))
saveas(gca, sprintf('%s/Phase_sinus_wave', fig_folder), 'epsc')

%% %% % % % %% % % %% % Fourier Transform % % % % % % % % % % % % %

Fs = 1e3; % sample rate. 1 Khz per second
T = 1/Fs; % Time steps
dt = 0:T:1; % duration 1 sec 

% signal frequencies
f1 = 10; 
f2 = 20;
f3 = 30;

s = 0.5*sin(2*pi*f1*dt') + 0.7*sin(2*pi*f2*dt') + sin(2*pi*f3*dt');
noisy_s = s + 2 .* normrnd(0,1, length(s),1);  % normally distributed noise

%
L = length(s);
ff = fft(s);  % Fourier transform of the signal
abs_norm_ff = abs(ff); %/max(abs(ff)); % normalisation
abs_norm_ff = abs_norm_ff/(length(ff)/2);

fff = ff(1:ceil(L/2))/L/2;
ff_noisy = fft(noisy_s); % fourier transform of noisy signal
abs_norm_ff_noisy = abs(ff_noisy)/(length(ff_noisy)/2); %/max(abs(ff_noisy));

figure
subplot(2,2,1)
plot(dt, s, 'LineWidth',2)
title('(a)', 'Interpreter','latex', 'FontSize',30)
xlabel('$t$ (sec)','Interpreter','latex', 'FontSize',20)
ylabel('$x(t)$', 'Interpreter','latex', 'FontSize',20)
ax=gca;
ax.FontSize = 30; 
yticks([-2, 0, 2])

subplot(2,2,3)
plot(dt, noisy_s, 'LineWidth',2)
title('(c)', 'Interpreter','latex', 'FontSize',30)
xlabel('$t$ (sec)','Interpreter','latex', 'FontSize',20)
ylabel('$x_{s}(t)$', 'Interpreter','latex', 'FontSize',20)
ax=gca;
ax.FontSize = 30; 
yticks([-10, 0, 10])
xticks()
% Fast Fourier Transformation on both signal and noisy signal. 

subplot(2,2,2)
stem((0:ceil(L/2)-1)*T*Fs, abs_norm_ff(1:ceil(end/2)), 'LineWidth',2);
% plot(f1,abs(ff))
title('(b)', 'Interpreter','latex', 'FontSize',30)
ylabel('$A(f)$', 'Interpreter','latex', 'FontSize',20)
xlabel('$f$ (Hz)', 'Interpreter','latex', 'FontSize',20)
ax=gca;
ax.FontSize = 30; 
xlim([0,50])

subplot(2,2,4)
stem((0:ceil(L/2)-1)*T*Fs, abs_norm_ff_noisy(1:ceil(end/2)), 'LineWidth',2);
% plot(f2,abs_norm_ff_noisy)
title('(d)', 'Interpreter','latex', 'FontSize',30)
ylabel('$A(f)$', 'Interpreter','latex', 'FontSize',20)
xlabel('$f$ (Hz)', 'Interpreter','latex', 'FontSize',20)
xlim([0,50])
ax=gca;
ax.FontSize = 30; 

%%

saveas(gca, sprintf('%s/FreqTimeDomainSine.fig', fig_folder))
saveas(gca, sprintf('%s/FreqTimeDomainSine', fig_folder), 'epsc')

%%
Adjacency= [1];
alpha=0;
x0 = rand(1,1);
r1 =3.5; % for regular trajectory
r2 = 4; % for fully-chaotic trajectory
niter= 1e4;

[Orbits_r1, Lyapunov_r1] = Log_generator(Adjacency, alpha, x0, r1, niter, 5e2, 1e-3, 1e4, 0);
[Orbits_r2, Lyapunov_r2] = Log_generator(Adjacency, alpha, x0, r2, niter, 5e2, 1e-3, 1e4, 0);
%
Fs = 100; % sample rate Khz
T = 1/Fs; % Time steps
dt = 0:T:niter-1; % Signal time duration

L1 = length(Orbits_r1); 
ff1_r1 = fft(Orbits_r1);
abs_norm_ff1_r1 = abs(ff1_r1)/(L1/2);

ff1_r2 = fft(Orbits_r2);
abs_norm_ff1_r2 = abs(ff1_r2)/(L1/2);
%%
figure
subplot(2,2,1)
plot(dt(1:50), Orbits_r1(1:50), 'LineWidth',2)
title('(a)', 'Interpreter','latex', 'FontSize',30)
ylim([0,1])
ax = gca;
ax.FontSize =30; 
ax.YTick = [0, 1];
ax.XTick = [0.2, 0.5];
ylabel("$x(t)$",'Interpreter','latex', 'FontSize',30 )
xlabel("$t$ (sec)",'Interpreter','latex', 'FontSize',30)

subplot(2,2,2)
plot((2:ceil(L1/2))*T, abs_norm_ff1_r1(2:ceil(end/2)), 'LineWidth',2);
title('(b)', 'Interpreter','latex', 'FontSize',30)
ylabel('$A(f)$', 'FontSize',30, 'Interpreter','latex')
xlabel('$f$ (Hz)', 'FontSize',30, 'Interpreter','latex')
ax = gca;
ax.FontSize =30;
ax.XTick = [20,50];
ax.YTick = [0, 0.06]; 

subplot(2,2,3)
plot(dt(1:50), Orbits_r2(1:50), 'LineWidth',2)
title('(c)', 'Interpreter','latex', 'FontSize',30)
ax = gca;
ax.FontSize =30;
% ax.YTick = [0, 1];
ax.XTick = [0.2, 0.5];
ax.YTick = [0, 1];
xlabel("$t$ (sec)",'Interpreter','latex', 'FontSize',30)
ylabel("$x(t)$",'Interpreter','latex', 'FontSize',30 )

subplot(2,2,4)
plot((2:ceil(L1/2))*T, abs_norm_ff1_r2(2:ceil(end/2)), 'LineWidth',2);
title('(d)', 'Interpreter','latex', 'FontSize',30)
ylabel('$A(f)$', 'FontSize',30, 'Interpreter','latex')
xlabel('$f$ (Hz)', 'FontSize',30, 'Interpreter','latex')
ax = gca;
ax.FontSize =30; ylim([0, 0.06]);
ax.YTick = [0, 0.06];
ax.XTick = [20,50];

%%
% 
saveas(gca, sprintf('%s/Fourier_Chaotic_Regular.fig', fig_folder))
saveas(gca, sprintf('%s/Fourier_Chaotic_Regular', fig_folder), 'epsc')

%% changing frequency in time, checking by Gabor transform

Fs = 1e2; % sample rate. 100 Hz per second
T = 1/Fs; % Time steps
dt = 0:T:Fs-T; % Signal time duration


% signal frequencies
f1 = 0.5.*dt; 

s = 0.5*sin(2*pi*f1.*dt);
L = length(s);
ff = fft(s);
abs_norm_ff = abs(ff)/(L/2);   % max(abs(ff));

%fff = ff(1:ceil(L/2))/L/2;
% ff_noisy = fft(noisy_s);
% abs_norm_ff_noisy = abs(ff_noisy)/max(abs(ff_noisy));

phs = hilbert_phases(s);
% phs = phs(1:100:end);
inst_freq = zeros(length(phs)-1, 1);
% %%
% L=length(s);
% C = fft(s);
% pow = (abs(C).^2)/L^2;
% 
% figure
% plot((0:ceil(L/2)-1)*T, pow2db(pow(1:ceil(end/2))/2).^2); 
% 
% %
% figure
% periodogram(s,rectwin(L), L, Fs, 'power')

%
for i=1:length(inst_freq)
    inst_freq(i) = (phs(i+1)-phs(i))/(2*pi*(T)); 
end

%%
figure
subplot(2,2,1)
plot(dt(1:1e3), s(1:1e3), 'LineWidth',2)
title('(a)', 'Interpreter','latex', 'FontSize',30)
xlabel('$t$ (sec)', 'FontSize',30, 'Interpreter','latex');
ylabel('$x_{t}$', 'FontSize',30, 'Interpreter','latex');
ax= gca; 
ax.FontSize = 30; 

subplot(2,2,2)
plot((0:ceil(L/2)-1)*T, abs_norm_ff(1:ceil(end/2)), 'LineWidth',2);
title('(b)', 'Interpreter','latex', 'FontSize',30)
xlabel('$f$ (Hz)', 'FontSize',30, 'Interpreter','latex');
ylabel('$A_{t}$', 'Interpreter','latex', 'FontSize',30)
ylim([0 1]); xlim([0,50])
ax= gca; 
ax.FontSize = 30; 
yticks([0, 0.5, 1]);

subplot(2,2,3)
spectrogram(s, 100, 16, 1e3, Fs,'yaxis'); % kaiser(128,18),120,128,1E3,'yaxis');

% spectrogram(s, 100, 80, 100, Fs,'yaxis'); % kaiser(128,18),120,128,1E3,'yaxis');
title('(c)', 'Interpreter','latex', 'FontSize',30)
xlabel('$t$ (min)', 'FontSize',30, 'Interpreter','latex')
% colormap jet
ax= gca; 
ax.FontSize = 30; 
ylabel('$f$ (Hz)', 'FontSize',30, 'interpreter', 'latex');
c1= colorbar;
% ylabel(c1,'Power (db)','FontSize',20,'Rotation',270)
c1.Limits = [-150,100];
c1.Ticks = [-150, 0, 100]; 

subplot(2,2,4)
plot(dt(2:end), inst_freq(1:end))
title('(d)', 'FontSize',30, 'Interpreter','latex')
ax= gca; 
ax.FontSize = 30; 
mintosecs = [0.5*60, 1*60, 1.5*60]; 
ax.XTick = mintosecs;
ax.XTickLabel = {'0.5', '1', '1.5'};
xlabel('$t$ (min)', 'FontSize',30, 'interpreter', 'latex')
ylabel('$f$ (Hz)', 'FontSize',30, 'interpreter', 'latex');
ylim([min(inst_freq), max(inst_freq)])
yticks([0, 20, 40])
%%

saveas(gca, sprintf('%s/Spectrogram_chirps.fig', fig_folder))
saveas(gca, sprintf('%s/Spectrogram_chirps', fig_folder), 'epsc')

%% Visualise HT for double sin wave

Fs = 1e3; % sample rate. 1 Khz per second
T = 1/Fs; % Time steps
dt = 0:T:Fs-T; % Signal time duration

f1 = 10; 
f2 = 20;
s = 0.5*sin(2*pi*f1*dt') + 0.5*sin(2*pi*f2*dt'); 

z = hilbert(s);
phs = unwrap(angle(z)); 
phs_range = linspace(0,2*pi,1e3); 
%% 
figure
subplot(2,2,[1:2]) % signal and envelope                        
plot(dt(1:1e3),s(1:1e3), 'LineWidth',2)
hold on
plot(dt(1:1e3), abs(z(1:1e3)), 'LineWidth',2, 'Color','#de7b12')
plot(dt(1:1e3), -abs(z(1:1e3)), 'LineWidth',2, 'Color','#de7b12')
ylim([-1.5,1.5])
ax = gca; 
ax.FontSize = 30; 
ax.XTick = [0, 0.5, 1];
xlabel('$t$ (sec)', 'Interpreter','latex', 'FontSize',30)
title('(a)', 'FontSize',30, 'Interpreter','latex')
legend('$x(t)$', '$\pm A(t)$', 'Interpreter','latex', 'FontSize',30)

subplot(2,2,3) % instantaneous phases vs. real phases
plot(dt(1:1e3), mod(phs(1:1e3),2*pi), 'LineWidth',2)
ax = gca; 
ax.FontSize = 30; 
ax.XTick = [0, 0.5, 1];
xlabel('$t$ (sec)', 'Interpreter','latex', 'FontSize',30)
ylabel('$\theta(t)$', 'Interpreter','latex', 'FontSize',30)
ax.YTick = [0, pi,  2*pi];
ax.YTickLabel = {'0', '\pi', '2\pi'};
ylim([0,2*pi])
title('(b)', 'FontSize',30, 'Interpreter','latex')

subplot(2,2,4)  % complex plane
% plot(cos(phs_range), sin(phs_range), '.', 'LineWidth',2)
pbaspect([1 1 1])
hold on 
plot(real(z(1:1e3)), imag(z(1:1e3)), 'LineWidth',2)
grid on
xlim([-2 2]); ylim([-2,2]); xline(0, 'LineWidth',1.5), yline(0, 'LineWidth',1.5);
hold on 
idx = find(round(phs,1) == round(pi/6,1));
plot([0, real(z(idx(1)))], [0, imag(z(idx(1)))], '.-','MarkerSize',16, 'LineWidth',2, 'color', "#e6742e"); 
text(0.3,0.1, '\pi/6', 'FontSize',20)
ax = gca; 
ax.FontSize = 30; 
title('(c)', 'FontSize',30, 'Interpreter','latex')
yticks([-2, 0, 2])
xlabel('$z_{R}(t)$', 'FontSize',30, 'Interpreter','latex')
ylabel('$z_{I}(t)$', 'FontSize',30, 'Interpreter','latex')

%%
saveas(gca, sprintf('%s/doublesinusoid.fig', fig_folder))
saveas(gca, sprintf('%s/doublesinusoid', fig_folder), 'epsc')













