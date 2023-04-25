function Get_fcut_quantiles(Occ_filename, Cuts_filename)
%% 16 April 2023. Script to obtain median data for different regions and Rpb3 deciles
% Includes Cuts 1 nt internal i.e. at G as well as A in GATC
% Divides by occupancy at G not A

% Site data are in the following .mat file
load('sacCer3_GATC_Site_Data.mat','ChrNoList','SiteList','Good_LeftSites','Good_RightSites',...
    'BadBothSites','ARS_Site','NDR_Site','ORF_Site','Rpb3_Decile',...
    'TEL_Site','Ty_Site','tRNA_Site');

%% Load Cuts and Occupancy data
load(Occ_filename, 'Occ');
load(Cuts_filename, 'LeftCut', 'RightCut');

%% Get fcut values for 2 nt at left and right sides of each GATC site. Use Occ at G and C
% Site coord is G in GATC. Left coord is the T; Right coord is the A (counter-intuitive!).
Right_Cut = zeros(35830,1);
RightOcc = zeros(35830,1);
Left_Cut = zeros(35830,1);
LeftOcc = zeros(35830,1);
Left_fcut = zeros(35830,1);
Right_fcut = zeros(35830,1);
for a = 1:35830
    Right_Cut(a) = RightCut{ChrNoList(a)}(SiteList(a)) + RightCut{ChrNoList(a)}(SiteList(a) +1);
    RightOcc(a) = Occ{ChrNoList(a)}(SiteList(a));
    Right_fcut(a) = Right_Cut(a)/RightOcc(a);
    Left_Cut(a) = LeftCut{ChrNoList(a)}(SiteList(a) +2) + LeftCut{ChrNoList(a)}(SiteList(a) +3);
    LeftOcc(a) = Occ{ChrNoList(a)}(SiteList(a) +3);
    Left_fcut(a) = Left_Cut(a)/LeftOcc(a);
end

% For sites < 200 bp apart on one side, left = right.
% For "bad" sites (< 200 bp away on both sides), = NaN.
for a = 1:35830
    if Good_LeftSites(a) == 0
     Left_fcut(a) = Right_fcut(a);
    end
    if Good_RightSites(a) == 0
     Right_fcut(a) = Left_fcut(a);
    end    
    if BadBothSites(a) == 1
     Left_fcut(a) = NaN;
     Right_fcut(a) = NaN;
    end 
end

%% Get median fcut for all good sites, using both left and right values
fcut = cat(1,Left_fcut, Right_fcut);
% Compute quantiles: 5%, 10%, 15%, etc.
fcut_Quantiles = quantile(fcut, 0.05:0.05:0.95)';

%% Get median fcut for good ORF, NDR, ARS, TEL, Ty sites
ORF_Left = zeros(35830,1);
ORF_Right = zeros(35830,1);
NDR_Left = zeros(35830,1);
NDR_Right = zeros(35830,1);
ARS_Left = zeros(35830,1);
ARS_Right = zeros(35830,1);
TEL_Left = zeros(35830,1);
TEL_Right = zeros(35830,1);
Ty_Left = zeros(35830,1);
Ty_Right = zeros(35830,1);
tRNA_Left = zeros(35830,1);
tRNA_Right = zeros(35830,1);

for a = 1:35830
    if ORF_Site(a) == 1
        ORF_Left(a) = Left_fcut(a);
        ORF_Right(a) = Right_fcut(a);
    elseif ORF_Site(a) == 0
        ORF_Left(a) = NaN;
        ORF_Right(a) = NaN;
    end
    if NDR_Site(a) == 1
       NDR_Left(a) = Left_fcut(a);
       NDR_Right(a) = Right_fcut(a);
    elseif NDR_Site(a) == 0
        NDR_Left(a) = NaN;
        NDR_Right(a) = NaN;
    end
    if ARS_Site(a) == 1
       ARS_Left(a) = Left_fcut(a);
       ARS_Right(a) = Right_fcut(a);
    elseif ARS_Site(a) == 0
        ARS_Left(a) = NaN;
        ARS_Right(a) = NaN;
    end
    if TEL_Site(a) == 1
       TEL_Left(a) = Left_fcut(a);
       TEL_Right(a) = Right_fcut(a);
    elseif TEL_Site(a) == 0
        TEL_Left(a) = NaN;
        TEL_Right(a) = NaN;
    end
    if Ty_Site(a) == 1
       Ty_Left(a) = Left_fcut(a);
       Ty_Right(a) = Right_fcut(a);
    elseif Ty_Site(a) == 0
        Ty_Left(a) = NaN;
        Ty_Right(a) = NaN;
    end
    if tRNA_Site(a) == 1
       tRNA_Left(a) = Left_fcut(a);
       tRNA_Right(a) = Right_fcut(a);
    elseif tRNA_Site(a) == 0
        tRNA_Left(a) = NaN;
        tRNA_Right(a) = NaN;
    end
end

% Calculate median for each data set and quantiles (10th quantile is median)
ORF_fcut = cat(1, ORF_Left, ORF_Right);
ORF_fcut_Quantiles = quantile(ORF_fcut, 0.05:0.05:0.95)';

NDR_fcut = cat(1, NDR_Left, NDR_Right);
NDR_fcut_Quantiles = quantile(NDR_fcut, 0.05:0.05:0.95)';

ARS_fcut = cat(1, ARS_Left, ARS_Right);
ARS_fcut_Quantiles = quantile(ARS_fcut, 0.05:0.05:0.95)';

TEL_fcut = cat(1, TEL_Left, TEL_Right);
TEL_fcut_Quantiles = quantile(TEL_fcut, 0.05:0.05:0.95)';

Ty_fcut = cat(1, Ty_Left, Ty_Right);
Ty_fcut_Quantiles = quantile(Ty_fcut, 0.05:0.05:0.95)';

tRNA_fcut = cat(1, tRNA_Left, tRNA_Right);
tRNA_fcut_Quantiles = quantile(tRNA_fcut, 0.05:0.05:0.95)';

%% CEN16 data: single site #34652, only left cut is good 
CEN16_fcut = Left_fcut(34652);

%% Medians and means for ORF deciles
Dec_ORF_Left = NaN(35830,10);
Dec_ORF_Right = NaN(35830,10);
Decile_Medians = zeros(10,1);
Decile_Means = zeros(10,1);

for a = 1:35830
    for b = 1:10
        if Rpb3_Decile(a) == b
        Dec_ORF_Left(a,b) = ORF_Left(a);
        Dec_ORF_Right(a,b) = ORF_Right(a);
        List = cat(1, Dec_ORF_Left(:,b), Dec_ORF_Right(:,b));
        Decile_Medians(b,1) = median(List,'omitnan');
        Decile_Means(b,1) = mean(List,'omitnan');
        end
    end
end

%% Save Quantile data (including median)
Q_filename = strrep(Cuts_filename,'Cuts_','Quant_');

save(sprintf(Q_filename),'fcut_Quantiles','ORF_fcut_Quantiles','NDR_fcut_Quantiles',...
    'ARS_fcut_Quantiles','TEL_fcut_Quantiles',...
    'Ty_fcut_Quantiles','tRNA_fcut_Quantiles','CEN16_fcut');

P_filename = strrep(Cuts_filename,'Cuts_','Rpb3_Deciles_');
save(sprintf(P_filename),'Decile_Medians','Decile_Means');

