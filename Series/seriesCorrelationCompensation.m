function [] = seriesCorrelationCompensation()
% removes the correlations from quadratures in mat-data, who are already
% piezoshaped. 

%% Discover *.mat files
filestruct = dir('mat-data/*.mat');
files = {filestruct.name};

%% Iterate through data files
dispstat('','init','timestamp','keepthis',0);

for i = 1:length(files)
    %% Load data
    C = strsplit(files{i},'.');
    filename = C{1};
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    clear rho WF;
    dispstat(['Loading ',files{i},' ...'],'timestamp',0);
    load(['mat-data/',files{i}]);
    [X1] = correlationCompensationAfterPiezo(X1);
    [X2] = correlationCompensationAfterPiezo(X2);
    [X3] = correlationCompensationAfterPiezo(X3);
    newFilename =  ['mat-data/',filename,'-corrRemove-yes.mat'];
    save(newFilename,'X1','X2','X3','piezoSign');
   
     
end

end



