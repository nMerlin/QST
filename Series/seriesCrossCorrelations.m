function [] = seriesCrossCorrelations(measurements,varargin)
%SERIESCROSSCORRELATIONS Plot cross-correlations for the given 3Ch-datasets

%% Validate and parse input arguments
p = inputParser;
defaultPaths = '.';
addParameter(p,'Paths',defaultPaths,@isvector);
defaultPrepOpts = struct;
addParameter(p,'PrepOpts',defaultPrepOpts,@isstruct);
defaultLOs = measurements-1;
addParameter(p,'LOs',defaultLOs,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[los,paths,prepopts] = c{:};

%% Create Cross-Correlation plots
for iPath = 1:length(paths) % iterate through different root paths
    cd(paths(iPath));
    listing = dir('raw-data');
    files = strjoin({listing.name}); % string containing all filenames
    for iMeas = 1:size(measurements,1) % iterate over measurements
        fLO = regexp(files, ...
            ['\d',num2str(los(iMeas)),'.*?.raw'],'match','once');
        fMeas = regexp(files, ...
            ['\d',num2str(measurements(iMeas)),'.*?.raw'],'match','once');
        [X1,X2,X3,~,config] = prepare3ChData(fLO,fMeas,prepopts);
        [A12,A13,A23] = plotCrossCorrelation(X1,X2,X3);
        % Add additional information to the plot
        n1 = var(X1(:)); f1 = config.E725.Piezo1.Frequency_DBL;
        n2 = var(X2(:)); f2 = config.E725.Piezo2.Frequency_DBL;
        n3 = var(X3(:)); f3 = config.E725.Piezo3.Frequency_DBL;
        title({['Smoothed Cross-Correlations for f1=',num2str(f1), ...
            ', f2=',num2str(f2),', f3=',num2str(f3),'Hz']; ...
            ['Photon numbers n1=',num2str(n1,3),', n2=',num2str(n2,3), ...
            ', n3=',num2str(n3,3)]; ...
            ['Amplitudes A12=',num2str(A12,3),', A13=',num2str(A13,3), ...
            ', A23=',num2str(A23,3)]});
        fOut = [datestr(date,'yyyy-mm-dd'),'-CrossCorrs-', ...
            fMeas,'.png'];
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [0 0 29.7 21];
        fig.PaperPositionMode = 'manual';
        print(fOut,'-dpng');
    end
end

