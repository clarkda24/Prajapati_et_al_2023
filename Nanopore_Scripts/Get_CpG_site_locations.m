function Get_CpG_site_locations(MSssI_filename)
%% 18 April 2023. Create table of CpG sites with locations
% Output - .csv and .mat files for site data, plus quantile, decile and phasing data

% CpG site data:
load(MSssI_filename,'Chr_No','Start_Coord','Central_Coord','End_Coord','Read_No', ...
    'Motif_No','m5C_No','Sequence','SssI_fmeth');

%% Identify sites in ORFs and get Rpb3 values (Pol II ChIP-seq data)
load('Rpb3_over_transcripts.mat','Avg_IP_over_Input');
load('Annotations_sacCer3.mat','Chr','ORF','ORFStart','ORFEnd','Watson');
Total_CpG_Sites = numel(Read_No);
ORF_Site = zeros(Total_CpG_Sites,1);
ORF_name = strings(Total_CpG_Sites,1);
Rpb3_level = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for h = 1:5770
        if Watson(h) == 1
            if Chr(h) == Chr_No(a)
                if Central_Coord(a) > ORFStart(h) && Central_Coord(a) < ORFEnd(h)
                ORF_Site(a) = 1;
                ORF_name(a) = ORF(h);
                Rpb3_level(a) = Avg_IP_over_Input(h);
                end
            end
            elseif Watson(h) == 0
                if Chr(h) == Chr_No(a)
                if Central_Coord(a) > ORFEnd(h) && Central_Coord(a) < ORFStart(h)
                ORF_Site(a) = 1; 
                ORF_name(a) = ORF(h);
                Rpb3_level(a) = Avg_IP_over_Input(h);
                end
                end
        end
    end
end

%% Assign Rpb3 decile - sort Rpb3 values in descending order
Rpb3_Deciles = sort(Avg_IP_over_Input, 'descend');
% Get the decile boundary values (1 to 5770) in the sorted lists 
Dvalues = Rpb3_Deciles([577,1154,1731,2308,2885,3462,4039,4616,5193]);
% Go through the CpG site list and determine which decile the value belongs to
Rpb3_Decile = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    if Rpb3_level(a) >= Dvalues(1,1)
        Rpb3_Decile(a) = 1;
    elseif Rpb3_level(a) < Dvalues(1,1) && Rpb3_level(a) >= Dvalues(2,1)
        Rpb3_Decile(a) = 2;
    elseif Rpb3_level(a) < Dvalues(2,1) && Rpb3_level(a) >= Dvalues(3,1)
        Rpb3_Decile(a) = 3;
    elseif Rpb3_level(a) < Dvalues(3,1) && Rpb3_level(a) >= Dvalues(4,1)
        Rpb3_Decile(a) = 4; 
    elseif Rpb3_level(a) < Dvalues(4,1) && Rpb3_level(a) >= Dvalues(5,1)
        Rpb3_Decile(a) = 5; 
    elseif Rpb3_level(a) < Dvalues(5,1) && Rpb3_level(a) >= Dvalues(6,1)
        Rpb3_Decile(a) = 6; 
    elseif Rpb3_level(a) < Dvalues(6,1) && Rpb3_level(a) >= Dvalues(7,1)
        Rpb3_Decile(a) = 7;  
    elseif Rpb3_level(a) < Dvalues(7,1) && Rpb3_level(a) >= Dvalues(8,1)
        Rpb3_Decile(a) = 8;
    elseif Rpb3_level(a) < Dvalues(8,1) && Rpb3_level(a) >= Dvalues(9,1)
        Rpb3_Decile(a) = 9;
    elseif Rpb3_level(a) < Dvalues(9,1) && Rpb3_level(a) > 0
        Rpb3_Decile(a) = 10;  
    end
end

%% Identify sites in promoter NDRs
% Plus1 and Minus1 are coordinates referring to the nucleosome dyads
% So NDR coords are from Minus1 dyad +73 to Plus1 dyad -73.
load('Annotations_sacCer3.mat','Plus1','Minus1','Watson');
NDR_Site = zeros(Total_CpG_Sites,1);
NDRwidth_check = zeros(5770,1);
for a = 1:Total_CpG_Sites
    for h = 1:5770
        if Watson(h) == 1
            if Chr(h) == Chr_No(a)
                if Central_Coord(a)+1 > Minus1(h) +73 && Central_Coord(a)+1 < Plus1(h) -73
                NDR_Site(a) = 1;
                NDRwidth_check(h) = Plus1(h) - Minus1(h) -146;
                end
            end
            elseif Watson(h) == 0
                if Chr(h) == Chr_No(a)
                    if Central_Coord(a)-1 > Plus1(h) +73 && Central_Coord(a)-1 < Minus1(h) -73
                    NDR_Site(a) = 1;
                    NDRwidth_check(h) = Minus1(h) - Plus1(h) -146;
                    end
                end
        end
    end
end

%% Identify sites in centromeres
load('CEN_data.mat','CEN_Chr_No','CEN_Left','CEN_Right','CEN_Top');
CEN_Site = zeros(Total_CpG_Sites,1);
CEN_Watson = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for d = 1:16
            if CEN_Chr_No(d) == Chr_No(a)
                if Central_Coord(a) >= CEN_Left(d) && Central_Coord(a) <= CEN_Right(d)
                CEN_Site(a) = 1;
                CEN_Watson(a) = CEN_Top(d);
                end
            end
    end
end

%% Identify ARS sites  
load('ARS_elements.mat','ARS_Chr_No','ARS_start','ARS_end');
ARS_Site = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for m = 1:352
            if ARS_Chr_No(m) == Chr_No(a)
                if Central_Coord(a) > ARS_start(m) && Central_Coord(a) < ARS_end(m)
                ARS_Site(a) = 1;
                end
            end
    end
end

%% Identify sites in telomeric regions (32 chromosome ends)
load('Telomeres.mat','TEL_Chr_No','TEL_start','TEL_end');
TEL_Site = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for r = 1:32
            if TEL_Chr_No(r) == Chr_No(a)
                if Central_Coord(a) > TEL_start(r) && Central_Coord(a) < TEL_end(r)
                TEL_Site(a) = 1;
                end
            end
    end
end

%% Identify sites in Ty elements and LTRs
% LTRs in Ty elements are also listed separately. 
load('Ty_elements.mat','Ty_Chr_No','Ty_start','Ty_end');
Ty_Site = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for s = 1:432
            if Ty_Chr_No(s) == Chr_No(a)
                if Central_Coord(a) > Ty_start(s) && Central_Coord(a) < Ty_end(s)
                Ty_Site(a) = 1;
                end
            end
    end
end

%% Identify sites in tRNA genes: 275 tRNA genes
clear 'Chr','Watson';
load('Sc3_tRNAgenes.mat','Chr','tRNAstart','tRNAend','Watson');
tRNA_Site = zeros(Total_CpG_Sites,1);
for a = 1:Total_CpG_Sites
    for j = 1:275
        if Watson(j) == 1
            if Chr(j) == Chr_No(a)
                if Central_Coord(a) > tRNAstart(j) && Central_Coord(a) < tRNAend(j)
                tRNA_Site(a) = 1;
                end
            end
            elseif Watson(j) == 0
                if Chr(j) == Chr_No(a)
                    if Central_Coord(a) > tRNAend(j) && Central_Coord(a) < tRNAstart(j)
                    tRNA_Site(a) = 1;
                    end
                end
        end
    end
end
clear 'Chr','Watson';

%% Save data as mat file
M_filename = strrep(MSssI_filename,'MSssI_','CpG_Site_Locations_');

save(sprintf(M_filename),'Chr_No','Start_Coord','End_Coord','Central_Coord',...
    'Motif_No','Sequence','Read_No','m5C_No','SssI_fmeth','ORF_name', ...
    'Rpb3_level','Rpb3_Decile','ORF_Site','NDR_Site','CEN_Site','CEN_Watson',...
    'ARS_Site','TEL_Site','Ty_Site','tRNA_Site');

%% A few CG sites have very low coverage because they are missing in our strain
% Define "Good_Cover" as > 10% of median:
Good_Cover = ones(numel(Read_No),1);
Median_Cover = median(Read_No(:), 'omitnan');
for n = 1:numel(Read_No)
    if Read_No(n) < 0.1 * Median_Cover
        Good_Cover(n) = 0;
    end
end

%% Get quantile data: group the CpG sites according to genomic region 
Total_Sites = numel(Chr_No);

fmeth_ORF = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if ORF_Site(a) == 1 && Good_Cover(a) == 1 
        fmeth_ORF(a) = SssI_fmeth(a);
    end
end
ORF_fmeth_Quantiles = quantile(fmeth_ORF,(0.05:0.05:0.95))';

fmeth_NDR = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if NDR_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_NDR(a) = SssI_fmeth(a);
    end
end
NDR_fmeth_Quantiles = quantile(fmeth_NDR,(0.05:0.05:0.95))';

fmeth_CEN = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if CEN_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_CEN(a) = SssI_fmeth(a);
    end
end
CEN_fmeth_Quantiles = quantile(fmeth_CEN,(0.05:0.05:0.95))';

fmeth_ARS = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if ARS_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_ARS(a) = SssI_fmeth(a);
    end
end
ARS_fmeth_Quantiles = quantile(fmeth_ARS,(0.05:0.05:0.95))';

fmeth_TEL = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if TEL_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_TEL(a) = SssI_fmeth(a);
    end
end
TEL_fmeth_Quantiles = quantile(fmeth_TEL,(0.05:0.05:0.95))';

fmeth_Ty = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if Ty_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_Ty(a) = SssI_fmeth(a);
    end
end
Ty_fmeth_Quantiles = quantile(fmeth_Ty,(0.05:0.05:0.95))';

fmeth_tRNA = NaN(Total_Sites,1);
for a = 1:Total_Sites
    if tRNA_Site(a) == 1 && Good_Cover(a) == 1
        fmeth_tRNA(a) = SssI_fmeth(a);
    end
end
tRNA_fmeth_Quantiles = quantile(fmeth_tRNA,(0.05:0.05:0.95))';

% Save quantile data
Q_filename = strrep(MSssI_filename,'MSssI_','Quantiles_CpG_');
save(sprintf(Q_filename),'ORF_fmeth_Quantiles','NDR_fmeth_Quantiles','CEN_fmeth_Quantiles',...
    'ARS_fmeth_Quantiles','TEL_fmeth_Quantiles','Ty_fmeth_Quantiles','tRNA_fmeth_Quantiles');

%% Get PolII Rpb3 decile data: means and medians
fmeth_Deciles = NaN(Total_Sites,10);
for a = 1:Total_Sites
    for b = 1:10
    if Rpb3_Decile(a) == b && Good_Cover(a) == 1
        fmeth_Deciles(a,b) = SssI_fmeth(a);
    end
    end
end
% Get mean and median for each decile
Decile_means = zeros(10,1);
Decile_medians = zeros(10,1);
for b = 1:10
    Decile_means(b,1) = mean(fmeth_Deciles(:,b),'omitnan');
    Decile_medians(b,1) = median(fmeth_Deciles(:,b),'omitnan');
end
% Save decile data
D_filename = strrep(MSssI_filename,'MSssI_','Deciles_CpG_');
save(sprintf(D_filename),'Decile_medians','Decile_means');

%% Phasing profiles: Get table of CpG site coords. rel. to Plus1. 
% Range = -1000 to +1000.
% Output: 'Phase_Coords_CpG_" .mat file for phasing plot
% Assign all sites relative to local +1 nucleosome coordinate = 0
% Collect all values at the same local coordinate - use cell arrays
Chromo = cell(1,16);
SiteCoord = cell(1,16);
PhaseCoord = cell(1,16);
Plus1Coord = cell(1,16);
WatCrick = cell(1,16);
fmeth = cell(1,16);

load('Annotations_sacCer3.mat','Chr','Watson');
% Phase coordinate corresponds to C or to the central nucleotide in a CpG cluster 
for h = 1:5770
    for a = 1:numel(Chr_No)
        if Chr(h) == Chr_No(a)
          if Watson(h) == 1
            if Central_Coord(a) - Plus1(h) < 1000 && Central_Coord(a) - Plus1(h) > -1000
            Chromo{Chr(h)}(end+1) = Chr_No(a);
            SiteCoord{Chr(h)}(end+1) = Central_Coord(a);
            PhaseCoord{Chr(h)}(end+1) = Central_Coord(a) - Plus1(h);
            Plus1Coord{Chr(h)}(end+1) = Plus1(h);
            WatCrick{Chr(h)}(end+1) = Watson(h);
            fmeth{Chr(h)}(end+1) = SssI_fmeth(a);
            end       
            elseif Watson(h) == 0
            if Plus1(h) - Central_Coord(a) < 1000 && Plus1(h) - Central_Coord(a) > -1000
            Chromo{Chr(h)}(end+1) = Chr_No(a);
            SiteCoord{Chr(h)}(end+1) = Central_Coord(a);
            PhaseCoord{Chr(h)}(end+1) = Plus1(h) - Central_Coord(a);
            Plus1Coord{Chr(h)}(end+1) = Plus1(h);
            WatCrick{Chr(h)}(end+1) = Watson(h);
            fmeth{Chr(h)}(end+1) = SssI_fmeth(a);
            end
          end
        end
    end
end

%% Construct mean phasing array: Local Coord (x) v. mean fmeth (y)
% Collect all fmeth values at the same local coordinate (-1000 to +1000)
fmeth = cell2mat(fmeth)';
Phase_Coord = cell2mat(PhaseCoord)';

fmeth_values = cell(1,2001);
for f = 1:numel(fmeth)
    for b = 1:2001
         if b - 1001 == Phase_Coord(f)
            fmeth_values{b}(end+1) = fmeth(f);
         end
    end
end
% Calculate mean at each local coordinate
fmeth_phase = cell(1,2001);
for c = 1:2001
    fmeth_phase{c} = mean(fmeth_values{c}(:),'omitnan');
end
fmeth_phase = cell2mat(fmeth_phase);
Sm_fmeth_phase = smooth(fmeth_phase, 21)';
% Save unsmoothed and smoothed phase profiles
S_filename = strrep(MSssI_filename,'MSssI_','Phase_Profiles_CpG_');
save(sprintf(S_filename),'fmeth_phase','Sm_fmeth_phase');

