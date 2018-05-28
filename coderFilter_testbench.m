% Testbench for the coderFilter function.
clear;
close all

%% Check dependencies

dependencies.toolboxDependencyAnalysis({'coderFilter'})

%% Sim settings
fs = 96e3;
dur = 1;
Ns = fs*dur;
f = (0:Ns-1).*(fs/Ns);
t = (0:Ns-1)./fs;

%% Signal setup
% Impulse response
x = [1 zeros(1,Ns-1)];

%% Cascade filters

% Composite band-pass filter
yDiscBand = coderFilter(x, fs, [10 1000], [2 1], {'high', 'low'});

% Increasing order low pass filter
nCascades = 4;
for mm=1:nCascades
    ynthLow(mm,:) = coderFilter(x, fs, 1e2*ones(1,mm),...
        1*ones(1,mm),...
        repmat({'low'},1,mm));
end

%% Plot

figure(1);
clf;
subplot(211)
for mm=1:nCascades
    semilogx(f, 20*log10(abs(fft(ynthLow(mm,:)))))
    hold on;
end
xlim([0 fs/2]);
ylim([-100 15]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend('2nd Order','4th Order', '6th Order', '8th Order');

subplot(212);
semilogx(f, 20*log10(abs(fft(yDiscBand))))
xlim([0 fs/2]);
ylim([-100 15]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend('Composite bandpass');

print('./img/example_response.png','-dpng','-r512');

%% MATLAB File Exchange logo

figure(2);
clf;

h = gcf;
set(h,'Units','Centimeters');
set(h, 'Position', [20 0 5 4])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[pos(3), pos(4)])


for mm=1:nCascades
    semilogx(f, 20*log10(abs(fft(ynthLow(mm,:)))),'LineWidth',1.5)
    hold on;
end
xlim([20 fs/2]);
ylim([-100 15]);
box off

set(gca,'xtick',[])
set(gca,'xticklabel',[])

set(gca,'ytick',[])
set(gca,'yticklabel',[])

print('./img/logo.png','-dpng','-r512');

