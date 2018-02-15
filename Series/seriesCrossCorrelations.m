function [] = seriesCrossCorrelations(measurements,varargin)
%SERIESCROSSCORRELATIONS Plot cross-correlations for the given 3Ch-datasets

%% Validate and parse input arguments
p = inputParser;
defaultPaths = {'.'};
addParameter(p,'Paths',defaultPaths,@isvector);
defaultPrepOpts = struct;
addParameter(p,'PrepOpts',defaultPrepOpts,@isstruct);
defaultLOs = measurements;
for iLOs=1:length(defaultLOs)
    defaultLOs{iLOs} = defaultLOs{iLOs} - 1;
end
addParameter(p,'LOs',defaultLOs,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[los,paths,prepopts] = c{:};

%% Create Cross-Correlation plots
dispstat('','init');
nPaths = length(paths);
for iPath = 1:nPaths % iterate through different root paths
    cd(paths{iPath});
    listing = dir('raw-data');
    files = strjoin({listing.name}); % string containing all filenames
    vMeas = measurements{iPath};
    vLos = los{iPath};
    nMeas = length(vMeas);
    for iMeas = 1:nMeas % iterate over measurements
        fLO = regexp(files, ...
            ['\d?',num2str(vLos(iMeas)),'-.*?.raw'],'match','once');
        fMeas = regexp(files, ...
            ['\d?',num2str(vMeas(iMeas)),'-.*?.raw'],'match','once');
        [X1,X2,X3,~,config] = prepare3ChData(fLO,fMeas,prepopts, ...
            'CorrRemove','yes');
        [A12,A13,A23] = plotCrossCorrelation(X1,X2,X3);
        % Add additional information to the plot
        n1 = var(X1(:))-0.5;
        n2 = var(X2(:))-0.5;
        n3 = var(X3(:))-0.5;
        % Read piezo frequencies
        if config.E725.Piezo1.GeneratorOn_BOOL
            f1 = config.E725.Piezo1.Frequency_DBL;
        else
            f1 = 0;
        end
        if config.E725.Piezo2.GeneratorOn_BOOL
            f2 = config.E725.Piezo2.Frequency_DBL;
        else
            f2 = 0;
        end
        if config.E725.Piezo3.GeneratorOn_BOOL
            f3 = config.E725.Piezo3.Frequency_DBL;
        else
            f3 = 0;
        end
        title({['Smoothed Cross-Correlations for f1=',num2str(f1), ...
            ', f2=',num2str(f2),', f3=',num2str(f3),'Hz']; ...
            ['Photon numbers n1=',num2str(n1,3),', n2=',num2str(n2,3), ...
            ', n3=',num2str(n3,3)]; ...
            ['Amplitudes A12=',num2str(A12,3),', A13=',num2str(A13,3), ...
            ', A23=',num2str(A23,3)]});
        fOut = [datestr(date,'yyyy-mm-dd'),'-CrossCorrs-CorrRemove-', ...
            fMeas,'.png'];
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [0 0 29.7 21];
        fig.PaperPositionMode = 'manual';
        print(fOut,'-dpng');
        dispstat(['Path ',num2str(iPath),' of ',num2str(nPaths), ...
            ' and measurement ',num2str(iMeas),' of ', ...
            num2str(nMeas),'.'],'keepthis');
        % Simulation
        nges = n1 + n2 + n3;
        r1 = sqrt(n1/nges);
        r2 = sqrt(n2/(nges-n1));
        [nPulses,nPieces,nSegments] = size(X1);
        [X1s,X2s,X3s] = simQST(nPulses,nPieces,nSegments,'Channels',3, ...
            'R1',r1,'R2',r2,'NPhotons',nges,'PhasePeriods',2.1, ...
            'Frequencies',[f1 f2 f3]);
        [A12s,A13s,A23s] = plotCrossCorrelation(X1s,X2s,X3s);
        title({['Simulated Smoothed Cross-Correlations for f1=',num2str(f1), ...
            ', f2=',num2str(f2),', f3=',num2str(f3),'Hz']; ...
            ['Photon numbers n1=',num2str(n1,3),', n2=',num2str(n2,3), ...
            ', n3=',num2str(n3,3)]; ...
            ['Amplitudes A12=',num2str(A12s,3),', A13=',num2str(A13s,3), ...
            ', A23=',num2str(A23s,3)]});
        fOut = [datestr(date,'yyyy-mm-dd'), ...
            '-Simulated-CrossCorrs-CorrRemove-',fMeas,'.png'];
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [0 0 29.7 21];
        fig.PaperPositionMode = 'manual';
        print(fOut,'-dpng');
    end
end

