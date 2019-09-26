
% Matlab code used for comparing the elemental signatures deposited at the
% vertebral edge of the pelagic juveniles and adults hammerhead captured 
% in the Mexican Pacific in 2016.

% By C. Coiraton October 2018

% All data manipulations and multivariate statistical analyses of were 
% performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
% Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
% College of Marine Science, University of South Florida, St. Petersburg, FL, USA
% http://www.marine.usf.edu/user/djones/matlab/matlab.html 


% Load data:
load all_HH.mat HH spots;

% Set up filename & save:
fname = 'HH_analysis_GEO_EDGE.Oct18';
saver;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  VARIABLES:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HH = structure with the following fields:
%    .oto     = numeric tag identifying individual otoliths
%    .txt     = cell array of text for each analyte element
%    .iso     = corresponding isotope number
%    .txt_iso = combined element + isotope labels
%    .nRep    = # within-otolith replicate spot samples used to calclate averages
%    .ppm     = concentration of unknown, averaged across replicates                 (ppm)
%    .LOD     = limits of detection, averaged across replicates                      (ppm)
%    .ratio   = molar ratios to internal standard, averaged across replicates (mMole/Mole)
%    .SRM     = name of SRM used for external calibration
%    .adj     = type of adjustment applied to values below LOD             (zero,LOD,none)
%    .spike   = type of spike removal applied to the time series
%    .drift   = method of drift correction
%    .tol     = tolerance for linear interpolation
%    .SRM     = name of SRM used for external calibration
%    .gDate   = cell array of Gregorian date of acquisition
% 
% 
% spots = structure with the following fields:
%   .spot   = unique spot scan identifier (average of 2-3 replicates)
%   .type   = type of spot sample (CO, PB, BM, GB, ED)
%   .loc    = location of collection site
%   .num    = specimen number from this site
%   .yr     = collection year
%   .mo     = collection month
%   .age    = age in years (-1 = embryo)
%   .sex    = F = female, M = male, E = embryo
%   .embryo = embryo number (0 = adult)
%   .slide  = master slide #
%   .vert   = number of vertebrae thin section on master slide
%   .gDate  = date of LA-ICP-MS run
% 
%   .loc_Rx  = location code
%   .type_Rx = (CO=1, PB=2, BM=3, GB=4, ED=5)
%   .GB_Rx   = growth band code (0-14)
%   .ED_Rx   = edge code (0-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  GEOGRAPHIC REGIONS (n = 4):                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create groups:
GRP_1 = {'BCS'}; % Baja California Sur (pelagic)
GRP_2 = {'MZT' }; % Mazatlán (coastal)
GRP_3 = {'ic'}; % Clarion Island (pelagic)
GRP_4 = {'PM'}; % Puerto Madero (coastal)

% Create grouping vector:
GEO = NaN(numel(spots.loc),1);
GEO(ismember(spots.loc,GRP_1)) = 1;
GEO(ismember(spots.loc,GRP_2)) = 2;
GEO(ismember(spots.loc,GRP_3)) = 3;
GEO(ismember(spots.loc,GRP_4)) = 4;
clear GRP*

% Verify groupings:
unique(spots.loc(GEO ==1))'
unique(spots.loc(GEO ==2))'
unique(spots.loc(GEO ==3))'
unique(spots.loc(GEO ==4))'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  COASTAL vs PELAGIC (n = 2):                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create groups
GRP_1 = {'MZT' 'PM' }; % COASTAL
GRP_2 = {'BCS' 'ic'}; %  %PELAGIC 

% Create grouping vector:
GEO = NaN(numel(spots.loc),1);
GEO(ismember(spots.loc,GRP_1)) = 1;
GEO(ismember(spots.loc,GRP_2)) = 2;
clear GRP*

% Verify groupings:
unique(spots.loc(GEO ==1))'
unique(spots.loc(GEO ==2))'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET INDIVIDUALS (ADULTS AND PELAGIC JUVENILES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get index to spot samples of ED (type_Rx = (CO=1, PB=2, BM=3, GB=4, ED=5)
idxED = (spots.ED_Rx==1);
sum(idxED)
idxA = (spots.age >= 5); %EXCLUDES EMBRYOS AS THEY ALSO HAVE EDGE SPOTS

% ALL INDIVIDUALS Get index to target locations:(LOC_RX: 2BCS 2OIC 13PM 12MZT 22ip) ALL
idxIC = ismember(spots.loc_Rx,[20]);

idx_MZT  = (spots.loc_Rx==12  & ismember(spots.num,[1 4 5])); %2016 only adults

idxPM  = (spots.loc_Rx==13  & ismember(spots.num,[6 7 8 9 10 11 12 13 14 15 ])); %PM 2016 only adults 

idxBCS = ismember(spots.loc_Rx,[2]);

idxLOC = (idxIC + idxPM + idxBCS + idxMZT);


% Get index to corresponding ED data:
idx = logical(idxED .* idxA .* idxLOC);
sum(idx)%48

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show group size:EDGE
sum(idx(GEO==1)) %16
sum(idx(GEO==2)) %3
sum(idx(GEO==3)) %19
sum(idx(GEO==4)) %10
sum(idx) %48

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data:
Y   = HH.ppm(idx,:);
LOD = HH.LOD(idx,:);

% Get percentage of each element below LOD:
bin  = Y<=LOD;
pLOD = (sum(bin)/size(bin,1));
pLOD'

% Get index to target elements (= not consistently below LOD):
tar = find(pLOD<0.1);

tar(ismember(tar,[4 5 12])) = []; % remove Ca43 (int.std, Ni60 (ICP cone made of Ni), P31(physiologically regulated)


% Confirm target elements are above LOD (should be all 0's):
Y(:,tar)<=LOD(:,tar)

% Get a list of the analytes:
HH.txt_iso(tar)'%list of the target analytes (<10% >LOD
LOD = HH.LOD(tar)' %mean LOD concentrations 

% 15 cell array
%     EDGE 
%     'Li7'
%     'Na23'
%     'Mg24'
%     'V51'
%     'Mn55'
%     'Fe57'
%     'Co59'
%     'Cu63'
%     'Zn64'
%     'Rb85'
%     'Sr88'
%     'Sn118'
%     'Ba137'
%     'Pb208'
%     'U238'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  PERMANOVA:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Test the null hypothesis of no difference in ED microchemistry signatures
% among the 4 geographic regions:

% Extract the relevant data:
Y   = HH.ppm(idx,tar);
grp = GEO(idx);

% Sort the data by group:
[~,key] = sort(grp); % get sort key
Y       = Y(key,:);  % sort Y
grp     = grp(key);  % sort grp

% Create Euclidean distance matrix from standardized data:
yDis = f_dis(f_stnd(Y),'euc');

%OPTIONAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----Variable Selection via RDA:-----
f_rdaStepwise(f_dummy(grp),Y,1000,1,0.05,0,0,HH.txt(tar))
% Remove elements that don't improve 'Cum R2adj':(USING RDA)
tar(ismember(tar,[13 17 ])) = []; % -> remove Rb, Cu  MAYBE REMOVE ZN CU AND MN
% %Update target elements:
Y = HH.ppm(idx,tar);
HH.txt_iso(tar)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: D.L. Jones wrote a new script (f_permanova and f_permanovaPW) in 
% May 2017 based off a paper published earlier that year by
% Anderson et al 2017 that internally 'transforms' residuals so that any
% combination of groups will fit this assumption. 

% PERMANOVA:
f_permanova(yDis,grp,1000,1);
% Perform pair-wise tests to see which groups differ:
f_permanovaPW(yDis,grp,1000,1);

% ==================================================
%                Modified PERMANOVA:            
% --------------------------------------------------
% F  = 2.9709   p    =  0.0010 (iter=1000) 
% ----------------------------------------------------------
% ----------------------------------------------------------
% Results of pair-wise comparisons between each factor level:
% ==========================================================
%           t:     p:     p_bon: p_ds:  p_holm:
% 1 vs. 2: 1.4942 0.0210 0.1260 0.1196 0.0420 
% 1 vs. 3: 1.8780 0.0020 0.0120 0.0119 0.0100 
% 1 vs. 4: 1.9615 0.0010 0.0060 0.0060 0.0060 
% 2 vs. 3: 1.5784 0.0270 0.1620 0.1515 0.0420 
% 2 vs. 4: 1.5200 0.0050 0.0300 0.0296 0.0150 
% 3 vs. 4: 1.8301 0.0030 0.0180 0.0179 0.0120 
% -------------------------------------------------------

%COASTAL VS PELAGIC
% ==================================================
%                Modified PERMANOVA:                
% --------------------------------------------------
%  F  = 3.2035    p    =  0.0010 (iter=1000) 
% --------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     CAP:                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Find optimal value of m:
f_capOptimal(f_stnd(Y),'euc',grp,1,1); %SPATIAL MEDIAN 
f_capOptimal(f_stnd(Y),'euc',grp,0,1); %CENTROID


% Perform CAP using optimal m:
cap = f_cap(f_stnd(Y),'euc',grp,[],1,1000,1,7,1); %SPATIAL MEDIAN ED  m=7 76.74 %    
cap = f_cap(f_stnd(Y),'euc',grp,[],1,1000,1,8,1); %COASTAL PELAGIC m=8 95.35%

cap = f_cap(f_stnd(Y),'euc',grp,[],0,1000,1,7,1);%=> CENTROID EDGE m=7 79.07 % 
cap = f_cap(f_stnd(Y),'euc',grp,[],0,1000,1,6,1);%%COASTAL PELAGIC m=6 93.02 % 


% CAP- EDGE SIGNATURE- SITES OF CAPTURE OF THE SHARKS- 

% ==================================================
%  CAP - Canonical Discriminant Analysis:
% --------------------------------------------------
% Trace Stat    = 1.6563  p =  0.00100 
% Greatest Root = 0.7590  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 7 
% Variability of yDis explained = 83.95 % 
% Canonical Correlations:
%   0.8712  0.7348  0.5979
% Squared Canonical Correlations (= delta^2):
%   0.7590  0.5399  0.3574
% ==================================================
% 
% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1             81.8 % 
%    2             66.7 % 
%    3             78.9 % 
%    4             80.0 % 
% 
% 
% Total Correct  = 79.07 % 
% Total Error    = 20.93 % 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      3      4      
%      1   81.8    0.0   18.2    0.0 
%      2    0.0   66.7    0.0   33.3 
%      3   21.1    0.0   78.9    0.0 
%      4    0.0    0.0   20.0   80.0 
% 
% ==================================================



% ==================================================
%  CAP - Canonical Discriminant Analysis:  COASTAL vs PELAGIC
% --------------------------------------------------
% Trace Stat    = 0.7565  p =  0.00100 
% Greatest Root = 0.7565  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 8 
% Variability of yDis explained = 87.95 % 
% Canonical Correlations:
%   0.8698
% Squared Canonical Correlations (= delta^2):
%   0.7565
% ==================================================
% 
% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1             84.6 % 
%    2            100.0 % 
% 
% 
% Total Correct  = 95.35 % 
% Total Error    = 4.65 % 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      
%      1   84.6   15.4 
%      2    0.0  100.0 
% 
% ==================================================


% Create CAP plot: what'a causing the differences, way to visualizer Permanova
% and find out more information about what the differences are and what's
% causing them

f_capPlot(cap,f_unique(num2str(grp)),[],Y,HH.txt(tar),0.03,'none',0,0,0,1);
f_capPlot(cap,f_unique(num2str(grp)),[],Y,HH.txt(tar),0.03,'none',0,0,1,1);%OVERLAyS VECTOR PLOT

% Save as PDF:
f_pdf('cap_4SITES');
f_pdf('cap_4SITES_vectors');

% Test the significance of the observed classification success rate:
f_chanceClass(grp,1-cap.loo_err.tot,1000,1);

% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1             9.88 % 
%    2             0.73 % 
%    3             9.88 % 
%    4             8.16 % 
% 
% Total Correct = 28.65 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 28.75 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------

% ==================================================
%            PROPORTIONAL CHANCE CRITERION: COASTAL VS PELAGIC
% --------------------------------------------------
% Group        Correct  
%    1            9.14 %  
%    2            48.67 %
% 
% Total Correct = 57.82 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 57.73 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                RANDOM FOREST:   (works with Matlab R2016b - not 2017                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Random Forests (RF; Breiman 2001) were employed in this study as non-linear 
%counterparts to the CAP classifiers described above. 
%RF were fitted using the same data set and evaluated similarly to estimate reclassification accuracies,
%construct confusion matrices using out-of-bag error rates (OOB-CV) and assess statistical
%significance using PCC. 
%This allowed a comparison of the results and performance of the two distinct methods
%for modeling vertebrae microchemistry data and the selection of the most accurate
%classifier to employ for habitat discrimination (Mercier et al. 2011)

% -----NOTES:----- (by David L. Jones)
% Used standardized variables (stnd = 1) to create a RF for DISCRIMINATION and
% assessing variable importance; use centered data (stnd = 2) to create a RF for
% CLASSIFYING unknowns.

rf_GEO = f_RFclass(Y,grp,1000,[],0,1,'stnd',1,HH.txt(tar));

% Test the significance of the observed classification success rate:
f_chanceClass(grp,0.6977,1000,1);

%  Create a Random Forest Canonical Discriminant Analysis plot:
f_RFvis(rf_GEO,0.95,0,0.05);

% Save as PDF:
f_pdf('rf_GEO');
f_pdf('rf_GEO_vectors');