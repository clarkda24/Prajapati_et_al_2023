function Get_fcut_phasing(Occ_filename, Cuts_filename)
%% 16 April 2023. Get fcut and plot phasing 
% 200 bp cut-off. Range -2000 to +2000. Plot: mean fcut from -1000 to +1000

% Data file
load('Phase_Coords.mat','Chromo_No','Site_Coord','Good_Left','Good_Right','Bad_Site',...
    'Phase_Coord');

%% Load Cuts and Occupancy data
load(Occ_filename, 'Occ');
load(Cuts_filename, 'LeftCut', 'RightCut');

%% Get fcut values for 2 nt at left and right sides of each GATC site. Use Occ at G and C
% Site coord is G in GATC. Left coord is the T; Right coord is the A (counter-intuitive!).
Right_Cut = zeros(66104,1);
RightOcc = zeros(66104,1);
Left_Cut = zeros(66104,1);
LeftOcc = zeros(66104,1);
Left_fcut = zeros(66104,1);
Right_fcut = zeros(66104,1);
for a = 1:66104
    Right_Cut(a) = RightCut{Chromo_No(a)}(Site_Coord(a)) + RightCut{Chromo_No(a)}(Site_Coord(a) +1);
    RightOcc(a) = Occ{Chromo_No(a)}(Site_Coord(a));
    Right_fcut(a) = Right_Cut(a)/RightOcc(a);
    Left_Cut(a) = LeftCut{Chromo_No(a)}(Site_Coord(a) +2) + LeftCut{Chromo_No(a)}(Site_Coord(a) +3);
    LeftOcc(a) = Occ{Chromo_No(a)}(Site_Coord(a) +3);
    Left_fcut(a) = Left_Cut(a)/LeftOcc(a);
end

% For sites < 200 bp apart on one side, left = right.
% For "bad" sites (< 200 bp away on both sides), = NaN.
for a = 1:66104
    if Good_Left(a) == 0
     Left_fcut(a) = Right_fcut(a);
    end
    if Good_Right(a) == 0
     Right_fcut(a) = Left_fcut(a);
    end    
    if Bad_Site(a) == 1
     Left_fcut(a) = NaN;
     Right_fcut(a) = NaN;
    end 
end

%% Construct mean phasing array: Local Coord (x) v. mean fcut (y)
% Collect all fcut values at the same local coordinate (-2000 to +2000)
fcut_values = cell(1,4001);
for a = 1:66104
    if Bad_Site(a) == 1
    end
        for b = 1:4001
                    if b - 2001 == Phase_Coord(a) 
                    fcut_values{b}(end+1) = Left_fcut(a);
                    fcut_values{b}(end+1) = Right_fcut(a);
                    end
        end
end
% Calculate mean at each local coordinate
fcut_phase = cell(1,4001);
for c = 1:4001
    fcut_phase{c} = mean(fcut_values{c}(:),'omitnan');
end
fcut_phase = cell2mat(fcut_phase);

% Plot mean fcut from -1000 to +1000
x = -1000:1000;
plot(x, fcut_phase(1000:3000));

%% Save mean fcut data
filename = strrep(Cuts_filename,'Cuts_','Mean_phase_fcut_');

save(sprintf(filename),'fcut_phase');

