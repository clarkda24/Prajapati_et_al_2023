%% 16 April 2023. Script to obtain phasing plot for 808H_2.
% Includes Cuts 1 nt internal i.e. at G as well as A in GATC
% Divides by occupancy at G not A
% Uses files produced by 'script_Get_phasing_data_808H_2.m'
% Three plots: unsmoothed, smoothed and smoothed +MNase-seq dyad data

Time = [0, 30, 60, 120, 240];

%% Create phase data array: unsmoothed and smoothed data
Phase_plot = zeros(5,4001);
Smoothed_plot = zeros(5,4001);

load('Mean_phase_fcut_808H_2_Dpn_0m_0_5000.mat');
Phase_plot(1,:) = fcut_phase(:);
Smoothed_plot(1,:) = smooth(fcut_phase(:),21);
clear 'fcut_phase';
load('Mean_phase_fcut_808H_2_Dpn_30m_0_5000.mat');
Phase_plot(2,:) = fcut_phase(:);
Smoothed_plot(2,:) = smooth(fcut_phase(:),21);
clear 'fcut_phase';
load('Mean_phase_fcut_808H_2_Dpn_60m_0_5000.mat');
Phase_plot(3,:) = fcut_phase(:);
Smoothed_plot(3,:) = smooth(fcut_phase(:),21);
clear 'fcut_phase';
load('Mean_phase_fcut_808H_2_Dpn_120m_0_5000.mat');
Phase_plot(4,:) = fcut_phase(:);
Smoothed_plot(4,:) = smooth(fcut_phase(:),21);
clear 'fcut_phase';
load('Mean_phase_fcut_808H_2_Dpn_240m_0_5000.mat');
Phase_plot(5,:) = fcut_phase(:);
Smoothed_plot(5,:) = smooth(fcut_phase(:),21);
clear 'fcut_phase';

%% Phasing plot -unsmoothed
figure('Position',[5,5,300,200]);
hold on
l(1) = plot(-1000:1000, Phase_plot(1,1000:3000));
l(2) = plot(-1000:1000, Phase_plot(2,1000:3000));
l(3) = plot(-1000:1000, Phase_plot(3,1000:3000));
l(4) = plot(-1000:1000, Phase_plot(4,1000:3000));
l(5) = plot(-1000:1000, Phase_plot(5,1000:3000));

legend(l, {'0 min','30','60','120', '240'}, 'location', 'EO','FontSize', 7)
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
xlabel('Position relative to +1 nucleosome (bp)', 'FontSize', 11)
title('808H_2 unsmoothed', 'interpreter', 'none', 'FontSize', 8)
set(gca, 'layer', 'top')
ylim([0, 1])
grid on
print(gcf, '-depsc', '-vector', 'Phasing_808H_2.eps');
hold off

%% Phasing plot -smoothed
figure('Position',[5,5,300,200]);
hold on
l(1) = plot(-1000:1000, Smoothed_plot(1,1000:3000));
l(2) = plot(-1000:1000, Smoothed_plot(2,1000:3000));
l(3) = plot(-1000:1000, Smoothed_plot(3,1000:3000));
l(4) = plot(-1000:1000, Smoothed_plot(4,1000:3000));
l(5) = plot(-1000:1000, Smoothed_plot(5,1000:3000));

legend(l, {'0 min','30','60','120', '240'}, 'location', 'EO','FontSize', 7)
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
xlabel('Position relative to +1 nucleosome (bp)', 'FontSize', 11)
title('808H_2 smoothed', 'interpreter', 'none', 'FontSize', 8)
set(gca, 'layer', 'top')
ylim([0, 1])
grid on
print(gcf, '-depsc', '-vector', 'Phasing_smoothed_808H_2.eps');
hold off

%% Phasing plot - smoothed, with MNase-seq dyads at 10%
load('Avg_dyad_density_WT_A_120_160_Ocampo_NAR_2016.mat', 'AvgDyads_Plus1');

figure('Position',[5,5,300,200]);
% For dyad area plot, grey fill:
area(-1000:1000, AvgDyads_Plus1/10, 'FaceColor', [0.9,0.9,0.9]);
hold on
l(1) = plot(-1000:1000, Smoothed_plot(1,1000:3000));
l(2) = plot(-1000:1000, Smoothed_plot(2,1000:3000));
l(3) = plot(-1000:1000, Smoothed_plot(3,1000:3000));
l(4) = plot(-1000:1000, Smoothed_plot(4,1000:3000));
l(5) = plot(-1000:1000, Smoothed_plot(5,1000:3000));

legend(l, {'0 min','30','60','120', '240'}, 'location', 'EO','FontSize', 7)
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
xlabel('Position relative to +1 nucleosome (bp)', 'FontSize', 11)
title('808H_2 smoothed', 'interpreter', 'none', 'FontSize', 8)
set(gca, 'layer', 'top')
ylim([0, 1])
grid on
print(gcf, '-depsc', '-vector', 'Phasing_smoothed_808H_2_dyads.eps');
hold off
