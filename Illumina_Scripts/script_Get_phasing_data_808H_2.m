%% 16 April 2023. Script to obtain data for the nucleosome phasing plot.
% Includes Cuts 1 nt internal i.e. at G as well as A in GATC
% Divides by occupancy at G not A
% Using function 'Get_fcut_phasing.m' 

Get_fcut_phasing('Occupancy_808H_2_Dpn_0m_0_5000.mat','Cuts_808H_2_Dpn_0m_0_5000.mat');

Get_fcut_phasing('Occupancy_808H_2_Dpn_30m_0_5000.mat','Cuts_808H_2_Dpn_30m_0_5000.mat');

Get_fcut_phasing('Occupancy_808H_2_Dpn_60m_0_5000.mat','Cuts_808H_2_Dpn_60m_0_5000.mat');

Get_fcut_phasing('Occupancy_808H_2_Dpn_120m_0_5000.mat','Cuts_808H_2_Dpn_120m_0_5000.mat');

Get_fcut_phasing('Occupancy_808H_2_Dpn_240m_0_5000.mat','Cuts_808H_2_Dpn_240m_0_5000.mat');

