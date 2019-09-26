% Examples for finding peaks or valleys in the elemental profile transects data of
% juveniles and adults hammerheads

% By C. Coiraton october 2018
% All data manipulations and multivariate statistical analyses of were 
% performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
% Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
% College of Marine Science, University of South Florida, St. Petersburg, FL, USA
% http://www.marine.usf.edu/user/djones/matlab/matlab.html 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inferring gabitat use changes and age related movements 
% based on Sr, Ba and Pb changes in vertebrae
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 34 - IC-LEW-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-14/170314.mat tra_34;
X = tra_34;

% Limit search to deeper and wider valleys:
[pks,loc] = f_peaks_PT(X,{'Pb208' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Ba137' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 31 - MZt-LEW-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-14/170314.mat tra_31;
X = tra_31;

% Limit search to deeper and wider valleys:
[pks,loc] = f_peaks_PT(X,{'Pb208' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Ba137' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 20 - PM-LEW-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-13/170313.mat tra_20;
X = tra_20;

% Limit search to deeper and wider valleys:
[pks,loc] = f_peaks_PT(X,{'Pb208' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Ba137' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 16 - BCS-LEW-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data:
load ./2017-03-11/170311.mat tra_16;
X = tra_16;

% Limit search to deeper and wider valleys:
[pks,loc] = f_peaks_PT(X,{'Pb208' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Ba137' 'Ca43'},1,0,11,1,3,2);
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%