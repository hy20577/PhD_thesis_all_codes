function  visualize_HT(Dat1, t, label)
 AS1= hilbert(Dat1);
% phase = hilbert_phases(Dat1);

%% %%% Visualisation of Hilbert Transform  (Double sinus wave)

% mp4 of phases and envelope 

fig = figure;
tiledlayout(2,1)
nexttile
plot(t, abs(Dat1))
h = animatedline(t(1), abs(AS1(1)), 'color', 'red', 'LineStyle','none');
% axis([0, 1.5, -pi, 5])
j = animatedline([t(1), t(1)], [0, abs(AS1(1))], 'color', 'red');
e = animatedline(t(1), abs(AS1(1)), 'color', 'red');
nexttile
g = animatedline([0 imag(AS1(1))],[0,real(AS1(1))], 'Marker','+', 'Color','red');
p = animatedline(imag(AS1(1)), real(AS1(1)), 'Color','b');
xline(0)
yline(0)
title('Analytical Signal', 'interpreter', 'latex', 'FontSize',20)
xlabel('Imaginary axis', 'interpreter', 'latex', 'FontSize',20)
ylabel('Real axis', 'interpreter', 'latex', 'FontSize',20)
pbaspect([1 1 1])
grid on
% axis([-4 4 -4 4])
set(gcf,'color','w'); % white background
% set(gcf, 'Position', get(0, 'Screensize'));
%%
obj = VideoWriter(sprintf('HT: %s', label), 'MPEG-4');
obj.Quality = 100;
obj.FrameRate = 40;
open(obj)

for k=2:length(Dat1)    
addpoints(h, t(k), abs(Dat1(k)));
clearpoints(j)
addpoints(j, [t(k), t(k)], [0, abs(AS1(k))])
addpoints(e, t(k), abs(AS1(k)))
clearpoints(g)
addpoints(g, [0 imag(AS1(k))],[0,real(AS1(k))])
addpoints(p, imag(AS1(k)), real(AS1(k)))
drawnow

frame= getframe(gcf);
writeVideo(obj, frame)

% frame for .png
    if k == ceil(length(Dat1)*0.6)
        fig_frame = getframe(gcf);
    end
end

obj.close();

figure
imshow(fig_frame.cdata)

saveas(gcf, sprintf('%s.fig', label))
saveas(gcf, sprintf('%s', label))
close






end