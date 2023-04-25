%% 16 April 2023. Median plots for various genomic regions and relative rate plot
% Includes Cuts 1 nt internal i.e. at G as well as A in GATC
% Divides by occupancy at G not A

% 5 time points
Time = [0, 30, 60, 120, 240];

% Collect the quantile data for each time point
ORF_Quant = zeros(5,19);
NDR_Quant = zeros(5,19);
ARS_Quant = zeros(5,19);
TEL_Quant = zeros(5,19);
Ty_Quant = zeros(5,19);
tRNA_Quant = zeros(5,19);
CEN16 = zeros(5,1);

load('Quant_808H_2_Dpn_0m_0_5000.mat');
ORF_Quant(1,:) = ORF_fcut_Quantiles(:);
NDR_Quant(1,:) = NDR_fcut_Quantiles(:);
ARS_Quant(1,:) = ARS_fcut_Quantiles(:);
TEL_Quant(1,:) = TEL_fcut_Quantiles(:);
Ty_Quant(1,:) = Ty_fcut_Quantiles(:);
tRNA_Quant(1,:) = tRNA_fcut_Quantiles(:);
CEN16(1,1) = CEN16_fcut;
clear 'fcut_Quantiles','ORF_fcut_Quantiles';'NDR_fcut_Quantiles';'ARS_fcut_Quantiles';...
    'TEL_fcut_Quantiles';'Ty_fcut_Quantiles';'tRNA_fcut_Quantiles';'CEN16_fcut';

load('Quant_808H_2_Dpn_30m_0_5000.mat');
ORF_Quant(2,:) = ORF_fcut_Quantiles(:);
NDR_Quant(2,:) = NDR_fcut_Quantiles(:);
ARS_Quant(2,:) = ARS_fcut_Quantiles(:);
TEL_Quant(2,:) = TEL_fcut_Quantiles(:);
Ty_Quant(2,:) = Ty_fcut_Quantiles(:);
tRNA_Quant(2,:) = tRNA_fcut_Quantiles(:);
CEN16(2,1) = CEN16_fcut;
clear 'fcut_Quantiles','ORF_fcut_Quantiles';'NDR_fcut_Quantiles';'ARS_fcut_Quantiles';...
    'TEL_fcut_Quantiles';'Ty_fcut_Quantiles';'tRNA_fcut_Quantiles';'CEN16_fcut';

load('Quant_808H_2_Dpn_60m_0_5000.mat');
ORF_Quant(3,:) = ORF_fcut_Quantiles(:);
NDR_Quant(3,:) = NDR_fcut_Quantiles(:);
ARS_Quant(3,:) = ARS_fcut_Quantiles(:);
TEL_Quant(3,:) = TEL_fcut_Quantiles(:);
Ty_Quant(3,:) = Ty_fcut_Quantiles(:);
tRNA_Quant(3,:) = tRNA_fcut_Quantiles(:);
CEN16(3,1) = CEN16_fcut;
clear 'fcut_Quantiles','ORF_fcut_Quantiles';'NDR_fcut_Quantiles';'ARS_fcut_Quantiles';...
    'TEL_fcut_Quantiles';'Ty_fcut_Quantiles';'tRNA_fcut_Quantiles';'CEN16_fcut';

load('Quant_808H_2_Dpn_120m_0_5000.mat');
ORF_Quant(4,:) = ORF_fcut_Quantiles(:);
NDR_Quant(4,:) = NDR_fcut_Quantiles(:);
ARS_Quant(4,:) = ARS_fcut_Quantiles(:);
TEL_Quant(4,:) = TEL_fcut_Quantiles(:);
Ty_Quant(4,:) = Ty_fcut_Quantiles(:);
tRNA_Quant(4,:) = tRNA_fcut_Quantiles(:);
CEN16(4,1) = CEN16_fcut;
clear 'fcut_Quantiles','ORF_fcut_Quantiles';'NDR_fcut_Quantiles';'ARS_fcut_Quantiles';...
    'TEL_fcut_Quantiles';'Ty_fcut_Quantiles';'tRNA_fcut_Quantiles';'CEN16_fcut';

load('Quant_808H_2_Dpn_240m_0_5000.mat');
ORF_Quant(5,:) = ORF_fcut_Quantiles(:);
NDR_Quant(5,:) = NDR_fcut_Quantiles(:);
ARS_Quant(5,:) = ARS_fcut_Quantiles(:);
TEL_Quant(5,:) = TEL_fcut_Quantiles(:);
Ty_Quant(5,:) = Ty_fcut_Quantiles(:);
tRNA_Quant(5,:) = tRNA_fcut_Quantiles(:);
CEN16(5,1) = CEN16_fcut;
clear 'fcut_Quantiles','ORF_fcut_Quantiles';'NDR_fcut_Quantiles';'ARS_fcut_Quantiles';...
    'TEL_fcut_Quantiles';'Ty_fcut_Quantiles';'tRNA_fcut_Quantiles';'CEN16_fcut';

%% Save median data (5 time points) for the 7 genomic regions
Medians = zeros(5,7);
Medians(:,1) = ORF_Quant(:,7);
Medians(:,2) = NDR_Quant(:,7);
Medians(:,3) = ARS_Quant(:,7);
Medians(:,4) = TEL_Quant(:,7);
Medians(:,5) = Ty_Quant(:,7);
Medians(:,6) = tRNA_Quant(:,7);
Medians(:,7) = CEN16(:,1);

% Save data as structure array
s1.Time = Time;
s1.Col_Names = {'ORF_ Sites', 'NDR_ Sites', 'ARS_ Sites',  ...
    'TEL_ Sites', 'Ty_ Sites', 'tRNA_ Sites','CEN16'};
s1.Medians = Medians;
save('808H_2_medians.mat','-struct','s1');

%% Plot for ORF good GATC sites (median = 10th quantile)
Quantiles = ORF_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_ORF', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6);
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_ORF.eps');
hold off
clear p

%% Plot for NDR good GATC sites (median = 10th quantile)
Quantiles = NDR_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_NDR', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6);  
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_NDR.eps');
hold off
clear p

%% Plot for ARS good GATC sites (median = 10th quantile)
Quantiles = ARS_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_ARS', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_ARS.eps');
hold off
clear p

%% Plot for TEL good GATC sites (median = 10th quantile)
Quantiles = TEL_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_TEL', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6);  
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_TEL.eps');
hold off
clear p

%% Plot for Ty good GATC sites (median = 10th quantile)
Quantiles = Ty_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_Ty', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_Ty.eps');
hold off
clear p

%% Plot for tRNA good GATC sites (median = 10th quantile)
Quantiles = tRNA_Quant';
% Colour settings
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Show percentile ranges
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_tRNA', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'SE','FontSize', 6);
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_tRNA.eps');
hold off
clear p

%% Plot for CEN16 - single site (left cut only)
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_CEN16', 'interpreter', 'none', 'FontSize', 8)
plot(x, CEN16, 'r', 'linewidth', 2); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_CEN16.eps');
hold off

%% Plots for Rpb3 deciles: Medians and Means
Decile_Med = zeros(5,10);
Decile_Av = zeros(5,10);

load('Rpb3_Deciles_808H_2_Dpn_0m_0_5000.mat');
Decile_Med(1,:) = Decile_Medians(:);
Decile_Av(1,:) = Decile_Means(:);
clear 'Decile_Medians'; 'Decile_Means';
load('Rpb3_Deciles_808H_2_Dpn_30m_0_5000.mat');
Decile_Med(2,:) = Decile_Medians(:);
Decile_Av(2,:) = Decile_Means(:);
clear 'Decile_Medians'; 'Decile_Means';
load('Rpb3_Deciles_808H_2_Dpn_60m_0_5000.mat');
Decile_Med(3,:) = Decile_Medians(:);
Decile_Av(3,:) = Decile_Means(:);
clear 'Decile_Medians'; 'Decile_Means';
load('Rpb3_Deciles_808H_2_Dpn_120m_0_5000.mat');
Decile_Med(4,:) = Decile_Medians(:);
Decile_Av(4,:) = Decile_Means(:);
clear 'Decile_Medians'; 'Decile_Means';
load('Rpb3_Deciles_808H_2_Dpn_240m_0_5000.mat');
Decile_Med(5,:) = Decile_Medians(:);
Decile_Av(5,:) = Decile_Means(:);
clear 'Decile_Medians'; 'Decile_Means';

% Plot for Rpb3 deciles MEDIANS
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Colour settings
newcolors = [1 0 0; 1 0 1; 1 1 0; 0.7 1 0; 0.3 1 0; 0 1 0; 0 1 1; 0 0.7 1; 0 0.3 1; 0 0 1];
colororder(newcolors)
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_Rpb3decile_Medians', 'interpreter', 'none', 'FontSize', 8)
l(1) = plot(x, Decile_Med(:,1), 'linewidth', 1);
l(2) = plot(x, Decile_Med(:,2), 'linewidth', 1);
l(3) = plot(x, Decile_Med(:,3), 'linewidth', 1);
l(4) = plot(x, Decile_Med(:,4), 'linewidth', 1);
l(5) = plot(x, Decile_Med(:,5), 'linewidth', 1);
l(6) = plot(x, Decile_Med(:,6), 'linewidth', 1);
l(7) = plot(x, Decile_Med(:,7), 'linewidth', 1);
l(8) = plot(x, Decile_Med(:,8), 'linewidth', 1);
l(9) = plot(x, Decile_Med(:,9), 'linewidth', 1);
l(10) = plot(x, Decile_Med(:,10), 'linewidth', 1);
legend({'Decile 1', 'Decile 2', 'Decile 3', 'Decile 4', 'Decile 5', 'Decile 6', 'Decile 7',...
    'Decile 8', 'Decile 9', 'Decile 10'}, 'FontSize', 6, 'location', 'SE'); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_808H_2_Rpb3deciles.eps');
hold off

% Plot for Rpb3 deciles MEANS
figure('Position', [5, 5, 225, 200])
x = Time;
hold all
% Colour settings
newcolors = [1 0 0; 1 0 1; 1 1 0; 0.7 1 0; 0.3 1 0; 0 1 0; 0 1 1; 0 0.7 1; 0 0.3 1; 0 0 1];
colororder(newcolors)
% Set X and Y limits
xlim([0 max(x)])
ylim([0 1])
% Set plot labels
xlabel('Time (min)', 'FontSize', 11)  
ylabel('Fraction methylated (cut by DpnI)', 'FontSize', 11)
title('808H_2_Rpb3decile_Means', 'interpreter', 'none', 'FontSize', 8)
l(1) = plot(x, Decile_Av(:,1), 'linewidth', 1);
l(2) = plot(x, Decile_Av(:,2), 'linewidth', 1);
l(3) = plot(x, Decile_Av(:,3), 'linewidth', 1);
l(4) = plot(x, Decile_Av(:,4), 'linewidth', 1);
l(5) = plot(x, Decile_Av(:,5), 'linewidth', 1);
l(6) = plot(x, Decile_Av(:,6), 'linewidth', 1);
l(7) = plot(x, Decile_Av(:,7), 'linewidth', 1);
l(8) = plot(x, Decile_Av(:,8), 'linewidth', 1);
l(9) = plot(x, Decile_Av(:,9), 'linewidth', 1);
l(10) = plot(x, Decile_Av(:,10), 'linewidth', 1);
legend({'Decile 1', 'Decile 2', 'Decile 3', 'Decile 4', 'Decile 5', 'Decile 6', 'Decile 7',...
    'Decile 8', 'Decile 9', 'Decile 10'}, 'FontSize', 6, 'location', 'SE');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Average_cut_ratios_808H_2_Rpb3deciles.eps');
hold off

%% Exponential plot for different regions.
% Plot ln(1-median fcut) for all regions v. time. Omit 0 time point
Log_Medians = log(1- Medians);

figure('Position', [5, 5, 225, 200])
hold all
xlim([0,250]);
xlabel('Time after SM addition (min)','FontSize', 11);
ylabel('ln(1 - fraction methylated)','FontSize', 11);
title('808H_2_Region_Log_Median', 'interpreter', 'none', 'FontSize', 8)
% Colour settings
NC = [1 0 0; 0 0 1; 1 0 1; 0.294 0.459 0.078; 0.85 0.325 0.098; 0 0.7 1; 0 0 0];
colororder(NC)

plot(Time(2:5), Log_Medians(2:5,1),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,2),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,3),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,4),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,5),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,6),'o', 'LineWidth', 1,'MarkerSize',2);
plot(Time(2:5), Log_Medians(2:5,7),'o', 'LineWidth', 1,'MarkerSize',2);
% Add least squares lines
lsline
set(gcf,'position',[5,5,225,200]);
legend('ORF','NDR','ARS','TEL','Ty','tRNA','CEN16','FontSize', 6','Location','SW');
hold off
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Log_median_808H_2.eps');
clear gcf

%% Calculate intercept and slope for log medians (y) v. time (x)
t = Time(2:5)';
tt = [ones(length(t),1) t];

% Calculate slopes of all 7 lines: Regression is X\Y 
Intercepts_Slopes = zeros(7,2);
for n = 1:7
    Intercept_Slope = tt\Log_Medians(2:5,n);
    Intercepts_Slopes(n,:) = Intercept_Slope;
end
Intercepts = Intercepts_Slopes(:,1);
Slopes = Intercepts_Slopes(:,2);

%% Calculate correlation coefficients for each plot predicted using y = mx +c
% R_Calc array: 4 time points, 7 sets of median data
R_Calc = zeros(4,7);
for m = 1:7
R_Calc(:,m) = tt(:,2) .* Slopes(m) + Intercepts(m);
end
% Calculate R-squared for each plot
% Rsq = 1 -[sum(y-yCalc),squared / sum(y- mean(y),squared)]
Rsq = zeros(7,1);
for p = 1:7
Rsq(p,1) = 1 - sum((Log_Medians(2:5,p) - R_Calc(:,p)).^2)/sum((Log_Medians(2:5,p) ...
    - mean(Log_Medians(2:5,p))).^2);
end

% Save the data in table form (.csv). 
Region = ["ORF","NDR","ARS","TEL","Ty","tRNA","CEN16"]';
T1 = table(Region, Slopes, Intercepts, Rsq, 'VariableNames',{'Region', ...
    'Slope','Intercept','R-squared'});
writetable(T1, 'Region_rates_plot_808H_2.csv');

