%% 18 April 2023. Median & Rel. Rate plots for various genomic regions: Rep 2.
% Includes rate plot (ln(1-regional_fmeth) v. time).
% 5 time points
Time = [0, 30, 60, 120, 240];

% Collect the quantile data for each time point
ORF_Quant = zeros(5,19);
NDR_Quant = zeros(5,19);
ARS_Quant = zeros(5,19);
TEL_Quant = zeros(5,19);
Ty_Quant = zeros(5,19);
tRNA_Quant = zeros(5,19);
CEN_Quant = zeros(5,19);

load('Quantiles_CpG_Rep2_0m.mat');
ORF_Quant(1,:) = ORF_fmeth_Quantiles(:);
NDR_Quant(1,:) = NDR_fmeth_Quantiles(:);
ARS_Quant(1,:) = ARS_fmeth_Quantiles(:);
TEL_Quant(1,:) = TEL_fmeth_Quantiles(:);
Ty_Quant(1,:) = Ty_fmeth_Quantiles(:);
tRNA_Quant(1,:) = tRNA_fmeth_Quantiles(:);
CEN_Quant(1,:) = CEN_fmeth_Quantiles(:);
clear 'ORF_fmeth_Quantiles';'NDR_fmeth_Quantiles';'ARS_fmeth_Quantiles';'TEL_fmeth_Quantiles'; ...
    'Ty_fmeth_Quantiles';'tRNA_fmeth_Quantiles';'CEN_fmeth_Quantiles';

load('Quantiles_CpG_Rep2_30m.mat');
ORF_Quant(2,:) = ORF_fmeth_Quantiles(:);
NDR_Quant(2,:) = NDR_fmeth_Quantiles(:);
ARS_Quant(2,:) = ARS_fmeth_Quantiles(:);
TEL_Quant(2,:) = TEL_fmeth_Quantiles(:);
Ty_Quant(2,:) = Ty_fmeth_Quantiles(:);
tRNA_Quant(2,:) = tRNA_fmeth_Quantiles(:);
CEN_Quant(2,:) = CEN_fmeth_Quantiles(:);
clear 'ORF_fmeth_Quantiles';'NDR_fmeth_Quantiles';'ARS_fmeth_Quantiles';'TEL_fmeth_Quantiles'; ...
    'Ty_fmeth_Quantiles';'tRNA_fmeth_Quantiles';'CEN_fmeth_Quantiles';

load('Quantiles_CpG_Rep2_60m.mat');
ORF_Quant(3,:) = ORF_fmeth_Quantiles(:);
NDR_Quant(3,:) = NDR_fmeth_Quantiles(:);
ARS_Quant(3,:) = ARS_fmeth_Quantiles(:);
TEL_Quant(3,:) = TEL_fmeth_Quantiles(:);
Ty_Quant(3,:) = Ty_fmeth_Quantiles(:);
tRNA_Quant(3,:) = tRNA_fmeth_Quantiles(:);
CEN_Quant(3,:) = CEN_fmeth_Quantiles(:);
clear 'ORF_fmeth_Quantiles';'NDR_fmeth_Quantiles';'ARS_fmeth_Quantiles';'TEL_fmeth_Quantiles'; ...
    'Ty_fmeth_Quantiles';'tRNA_fmeth_Quantiles';'CEN_fmeth_Quantiles';

load('Quantiles_CpG_Rep2_120m.mat');
ORF_Quant(4,:) = ORF_fmeth_Quantiles(:);
NDR_Quant(4,:) = NDR_fmeth_Quantiles(:);
ARS_Quant(4,:) = ARS_fmeth_Quantiles(:);
TEL_Quant(4,:) = TEL_fmeth_Quantiles(:);
Ty_Quant(4,:) = Ty_fmeth_Quantiles(:);
tRNA_Quant(4,:) = tRNA_fmeth_Quantiles(:);
CEN_Quant(4,:) = CEN_fmeth_Quantiles(:);
clear 'ORF_fmeth_Quantiles';'NDR_fmeth_Quantiles';'ARS_fmeth_Quantiles';'TEL_fmeth_Quantiles'; ...
    'Ty_fmeth_Quantiles';'tRNA_fmeth_Quantiles';'CEN_fmeth_Quantiles';

load('Quantiles_CpG_Rep2_240m.mat');
ORF_Quant(5,:) = ORF_fmeth_Quantiles(:);
NDR_Quant(5,:) = NDR_fmeth_Quantiles(:);
ARS_Quant(5,:) = ARS_fmeth_Quantiles(:);
TEL_Quant(5,:) = TEL_fmeth_Quantiles(:);
Ty_Quant(5,:) = Ty_fmeth_Quantiles(:);
tRNA_Quant(5,:) = tRNA_fmeth_Quantiles(:);
CEN_Quant(5,:) = CEN_fmeth_Quantiles(:);
clear 'ORF_fmeth_Quantiles';'NDR_fmeth_Quantiles';'ARS_fmeth_Quantiles';'TEL_fmeth_Quantiles'; ...
    'Ty_fmeth_Quantiles';'tRNA_fmeth_Quantiles';'CEN_fmeth_Quantiles';

%% Plot for ORF good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - ORF', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW','FontSize', 6);   
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_ORF.eps');
hold off
clear p

%% Plot for NDR good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - NDR', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW','FontSize', 6);   
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_NDR.eps');
hold off
clear p

%% Plot for ARS good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - ARS', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW', 'FontSize', 6); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_ARS.eps');
hold off
clear p

%% Plot for TEL good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - TEL', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW','FontSize', 6);  
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_TEL.eps');
hold off
clear p

%% Plot for Ty good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - Ty', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW','FontSize', 6);
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_Ty.eps');
hold off
clear p

%% Plot for tRNA good CpG sites (median = 10th quantile)
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - tRNA', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW', 'FontSize', 6);   
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_tRNA.eps');
hold off
clear p

%% Plot for centromeres - good CpG sites (median = 10th quantile)
Quantiles = CEN_Quant';
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation by M.SssI, Rep2 - CEN', 'interpreter', 'none', 'FontSize', 8)
plot(x, Quantiles(10,:), 'r', 'linewidth', 2);
h = legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'NW', 'FontSize', 6); 
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_cut_ratios_SssI_Rep2_CEN.eps');
hold off
clear p

%% Plots for Rpb3 deciles: Medians and Means
Decile_Med = zeros(5,10);
Decile_Av = zeros(5,10);

load('Deciles_CpG_Rep2_0m.mat');
Decile_Med(1,:) = Decile_medians(:);
Decile_Av(1,:) = Decile_means(:);
clear 'Decile_medians'; 'Decile_means';
load('Deciles_CpG_Rep2_30m.mat');
Decile_Med(2,:) = Decile_medians(:);
Decile_Av(2,:) = Decile_means(:);
clear 'Decile_medians'; 'Decile_means';
load('Deciles_CpG_Rep2_60m.mat');
Decile_Med(3,:) = Decile_medians(:);
Decile_Av(3,:) = Decile_means(:);
clear 'Decile_medians'; 'Decile_means';
load('Deciles_CpG_Rep2_120m.mat');
Decile_Med(4,:) = Decile_medians(:);
Decile_Av(4,:) = Decile_means(:);
clear 'Decile_medians'; 'Decile_means';
load('Deciles_CpG_Rep2_240m.mat');
Decile_Med(5,:) = Decile_medians(:);
Decile_Av(5,:) = Decile_means(:);
clear 'Decile_medians'; 'Decile_means';

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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation, Rep2: Rpb3 Decile Medians', 'interpreter', 'none', 'FontSize', 8)
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
    'Decile 8', 'Decile 9', 'Decile 10'}, 'FontSize', 10, 'location', 'NW', 'FontSize', 6);  
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Median_CpG_methylation_Rpb3_deciles_Rep2.eps');
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
xlabel('Time after SM addition (min)', 'FontSize', 11)  
ylabel('CpG fraction methylated', 'FontSize', 11)
title('CpG methylation, Rep2: Rpb3 Decile Means', 'interpreter', 'none', 'FontSize', 8)
l(1) = plot(x, Decile_Av(:,1), 'linewidth', 1);
l(2) = plot(x, Decile_Av(:,2), 'linewidth', 1);
l(3) = plot(x, Decile_Av(:,3), 'linewidth', 1);
l(4) = plot(x, Decile_Av(:,4), 'linewidth', 1);
l(5) = plot(x, Decile_Av(:,5), 'linewidth', 1);
l(6) = plot(x, Decile_Av(:,6), 'linewidth', 1);
l(7) = plot(x, Decile_Av(:,7), 'linewidth', 1);
l(8) = plot(x, Decile_Av(:,8), 'linewidth', 1);
l(9) = plot(x, Decile_Av(:,9), 'linewidth', 1);
l(10) = plot(x, Decile_Av(:,10), 'linewidth', 1);legend({'Decile 1', 'Decile 2', 'Decile 3', 'Decile 4', 'Decile 5', 'Decile 6', 'Decile 7',...
    'Decile 8', 'Decile 9', 'Decile 10'}, 'FontSize', 10, 'location', 'NW', 'FontSize', 6);  
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Mean_CpG_methylation_Rpb3_deciles_Rep2.eps');
hold off

%% Save median data (5 time points) for the 7 genomic regions
Medians = zeros(5,7);
Medians(:,1) = ORF_Quant(:,7);
Medians(:,2) = NDR_Quant(:,7);
Medians(:,3) = ARS_Quant(:,7);
Medians(:,4) = TEL_Quant(:,7);
Medians(:,5) = Ty_Quant(:,7);
Medians(:,6) = tRNA_Quant(:,7);
Medians(:,7) = CEN_Quant(:,7);

%% Plot and save exponential decay data for different regions.
% Plot ln(1-median fcut) for all regions v. time. Omit 0 time point
Log_Medians = log(1- Medians);

figure('Position', [5, 5, 225, 200])
hold all
xlim([0,250]);
ylim([-2.5,0.5]);
xlabel('Time after SM addition (min)','FontSize', 11);
ylabel('ln(1 - fraction methylated)','FontSize', 11);
title('Log Median MSssI Rep2', 'interpreter', 'none', 'FontSize', 8)
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
legend('ORF','NDR','ARS','TEL','Ty','tRNA','CEN','FontSize', 6','Location','SW');
hold off
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
print(gcf, '-depsc', '-vector', 'Log_medians_MSssI_Rep2.eps');
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
for m = 1:7
Rsq(m,1) = 1 - sum((Log_Medians(2:5,m) - R_Calc(:,m)).^2)/sum((Log_Medians(2:5,m) ...
    - mean(Log_Medians(2:5,m))).^2);
end

% Save the data in table form (.csv). 
Region = ["ORF","NDR","ARS","TEL","Ty","tRNA","CEN"]';
T1 = table(Region, Slopes, Intercepts, Rsq,'VariableNames',{'Region', ...
    'Slope','Intercept','R-squared'});
writetable(T1, 'Region_rates_plot_MSssI_Rep2.csv');

