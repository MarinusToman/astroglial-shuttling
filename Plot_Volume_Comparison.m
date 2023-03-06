close all;
clear all variables;

fileSml = 'Astro_Shuttling_Model_SynCovA_10Hz_23_Jan_2023.fig';
fileNml = 'Astro_Shuttling_Model_SynCovB_10Hz_23_Jan_2023.fig';
fileBig = 'Astro_Shuttling_Model_SynCovC_10Hz_23_Jan_2023.fig';

%% Small volume
figSmall = openfig(fileSml);% get figure handle
plotHandle = findobj(figSmall(3),'Type','axes'); %get third figure handle

plotLine = findobj(plotHandle(16),'Type','line');% get K_psc subplot
t = plotLine.XData;% time data
K_psc_small = plotLine.YData; % K_psc data
plotLine = findobj(plotHandle(10),'Type','line');% get K_ecs subplot
K_ecs_small = plotLine.YData; % K_ecs data

plotLine = findobj(plotHandle(18),'Type','line');% get Na_psc subplot
Na_psc_small = plotLine.YData; % Na_psc data
plotLine = findobj(plotHandle(12),'Type','line');% get Na_ecs subplot
Na_ecs_small = plotLine.YData; % Na_ecs data

plotLine = findobj(plotHandle(14),'Type','line');% get Ca_psc subplot
Ca_psc_small = plotLine.YData; % Ca_psc data
plotLine = findobj(plotHandle(8),'Type','line');% get Ca_ecs subplot
Ca_ecs_small = plotLine.YData; % Ca_ecs data

close all;

%% Normal volume
fig = openfig(fileNml);% get figure handle
plotHandle = findobj(fig(3),'Type','axes'); %get third figure handle

plotLine = findobj(plotHandle(16),'Type','line');% get K_psc subplot
K_psc = plotLine.YData; % K_psc data
plotLine = findobj(plotHandle(10),'Type','line');% get K_ecs subplot
K_ecs = plotLine.YData; % K_ecs data

plotLine = findobj(plotHandle(18),'Type','line');% get Na_psc subplot
Na_psc = plotLine.YData; % Na_psc data
plotLine = findobj(plotHandle(12),'Type','line');% get Na_ecs subplot
Na_ecs = plotLine.YData; % Na_ecs data

plotLine = findobj(plotHandle(14),'Type','line');% get Ca_psc subplot
Ca_psc = plotLine.YData; % Ca_psc data
plotLine = findobj(plotHandle(8),'Type','line');% get Ca_ecs subplot
Ca_ecs = plotLine.YData; % Ca_ecs data

close all;

%% Big volume
figBig = openfig(fileBig);% get figure handle
plotHandle = findobj(figBig(3),'Type','axes'); %get third figure handle

plotLine = findobj(plotHandle(16),'Type','line');% get K_psc subplot
K_psc_big = plotLine.YData; % K_psc data
plotLine = findobj(plotHandle(10),'Type','line');% get K_ecs subplot
K_ecs_big = plotLine.YData; % K_ecs data

plotLine = findobj(plotHandle(18),'Type','line');% get Na_psc subplot
Na_psc_big = plotLine.YData; % Na_psc data
plotLine = findobj(plotHandle(12),'Type','line');% get Na_ecs subplot
Na_ecs_big = plotLine.YData; % Na_ecs data

plotLine = findobj(plotHandle(14),'Type','line');% get Ca_psc subplot
Ca_psc_big = plotLine.YData; % Ca_psc data
plotLine = findobj(plotHandle(8),'Type','line');% get Ca_ecs subplot
Ca_ecs_big = plotLine.YData; % Ca_ecs data

%% Plot Data
close all;
fig = figure('Name','Concentrations in Different ECS Volumes');
fig.Units = 'normalized';
fig.OuterPosition = [0 0 1 1];

blCol = [0,0.4470,0.7410];
yCol = [0.9290 0.6940 0.1250];
oCol = [0.8500 0.3250 0.0980];

subplot(3,3,1)
plot(t,Na_psc_small,'Color',yCol);hold on;
plot(t,Na_psc,'Color',blCol);
plot(t,Na_psc_big,'Color',oCol);hold off;
title('[Na^+]_{PsC}');
ylabel('[Na^+]_{PsC} (mM)');
legend('A','B','C');

subplot(3,3,2)
plot(t,Na_ecs_small,'Color',yCol);hold on;
plot(t,Na_ecs,'Color',blCol);
plot(t,Na_ecs_big,'Color',oCol);hold off;
title('[Na^+]_{ECS}');
ylabel('[Na^+]_{ECS} (mM)');

subplot(3,3,4)
plot(t,K_psc_small,'Color',yCol);hold on;
plot(t,K_psc,'Color',blCol);
plot(t,K_psc_big,'Color',oCol);hold off;
title('[K^+]_{PsC}');
ylabel('[K^+]_{PsC} (mM)');

subplot(3,3,5)
plot(t,K_ecs_small,'Color',yCol);hold on;
plot(t,K_ecs,'Color',blCol);
plot(t,K_ecs_big,'Color',oCol);hold off;
title('[K^+]_{ECS}');
ylabel('[K^+]_{ECS} (mM)');

subplot(3,3,7)
plot(t,Ca_psc_small,'Color',yCol);hold on;
plot(t,Ca_psc,'Color',blCol);
plot(t,Ca_psc_big,'Color',oCol);hold off;
title('[Ca^{2+}]_{PsC}');
ylabel('[Ca^{2+}]_{PsC} (nM)');
xlabel('Time (s)');

subplot(3,3,8)
plot(t,Ca_ecs_small,'Color',yCol);hold on;
plot(t,Ca_ecs,'Color',blCol);
plot(t,Ca_ecs_big,'Color',oCol);hold off;
title('[Ca^{2+}]_{ECS}');
ylabel('[Ca^{2+}]_{ECS} (mM)');
xlabel('Time (s)');
% legend('A','B','C');

savefig(fig,'Compare_Synapse_Coverage');