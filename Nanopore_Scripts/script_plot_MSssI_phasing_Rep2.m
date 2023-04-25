%% 18 April 23. Phasing plots for M.SssI data Rep 2. 
% Mean fmeth from -1000 to +1000 (unsmoothed and smoothed)
% With MNase-seq dyad data
load('Phase_Profiles_CpG_y111_Con2.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_con = fmeth_phase;
Sm_fmeth_phase_con = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Phase_Profiles_CpG_Rep2_0m.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_0m = fmeth_phase;
Sm_fmeth_phase_0m = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Phase_Profiles_CpG_Rep2_30m.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_30m = fmeth_phase;
Sm_fmeth_phase_30m = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Phase_Profiles_CpG_Rep2_60m.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_60m = fmeth_phase;
Sm_fmeth_phase_60m = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Phase_Profiles_CpG_Rep2_120m.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_120m = fmeth_phase;
Sm_fmeth_phase_120m = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Phase_Profiles_CpG_Rep2_240m.mat','fmeth_phase','Sm_fmeth_phase');
fmeth_phase_240m = fmeth_phase;
Sm_fmeth_phase_240m = Sm_fmeth_phase;
clear fmeth_phase Sm_fmeth_phase

load('Avg_dyad_density_WT_A_120_160_Ocampo_NAR_2016.mat', 'AvgDyads_Plus1');

%% Plot unsmoothed profiles 
figure('Position', [5, 5, 300, 200])
% For dyad area plot, grey fill:
area(-1000:1000, AvgDyads_Plus1/10, 'FaceColor', [0.9,0.9,0.9]);
hold on
x = -1000:1000;
plot(x, fmeth_phase_con,'r');
ylim([0, 1])
plot(x, fmeth_phase_0m,'m');
plot(x, fmeth_phase_30m,'g');
plot(x, fmeth_phase_60m,'c');
plot(x, fmeth_phase_120m,'b');
plot(x, fmeth_phase_240m,'k');
legend({'MNase dyads','Control','0 min','30 min','60 min','120 min','240 min'},...
    'location', 'EO','FontSize', 7);
ylabel('Fraction methylated by M.SssI', 'FontSize', 11)
xlabel('Position relative to +1 nucleosome (bp)', 'FontSize', 11)
title('Phasing M.SssI Rep 2, unsmoothed', 'interpreter', 'none', 'FontSize', 8)
set(gca, 'layer', 'top')
grid on
print(gcf, '-depsc', '-vector', 'Phasing_M.SssI_Rep2_dyads_unsmoothed.eps');
hold off
clear gcf

%% Plot smoothed profiles 
figure('Position', [5, 5, 300, 200])
% For dyad area plot, grey fill:
area(-1000:1000, AvgDyads_Plus1/10, 'FaceColor', [0.9,0.9,0.9]);
hold on
x = -1000:1000;
plot(x, Sm_fmeth_phase_con,'r','LineWidth',1);
ylim([0, 1])
plot(x, Sm_fmeth_phase_0m,'m','LineWidth',1);
plot(x, Sm_fmeth_phase_30m,'g','LineWidth',1);
plot(x, Sm_fmeth_phase_60m,'c','LineWidth',1);
plot(x, Sm_fmeth_phase_120m,'b','LineWidth',1);
plot(x, Sm_fmeth_phase_240m,'k','LineWidth',1);
legend({'MNase dyads','Control','0 min','30 min','60 min','120 min','240 min'},...
    'location', 'EO','FontSize', 7);
ylabel('Fraction methylated by M.SssI', 'FontSize', 11)
xlabel('Position relative to +1 nucleosome (bp)', 'FontSize', 11)
title('Phasing M.SssI Rep 2, smoothed', 'interpreter', 'none', 'FontSize', 8)
set(gca, 'layer', 'top')
grid on
print(gcf, '-depsc', '-vector', 'Phasing_M.SssI_Rep2_dyads_smoothed.eps');
hold off
clear gcf
