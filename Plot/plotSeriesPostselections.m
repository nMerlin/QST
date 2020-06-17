function plotSeriesPostselections(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultType = 'Amplitude';
addParameter(p,'Type',defaultType,@isstr);
defaultXUnit = 'fs';
addParameter(p,'XUnit',defaultXUnit,@isstr);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,typestr,varyAPS,xUnit] = c{:};

% Constants
figurepath = 'figures-fig/';

%% Gather data
[delay,Yr,Yt,discAmpl,discMeanVar,discN,g2vals,g2std] = deal([]);
sigmas = zeros(length(listOfParams),1);
sigmaConf = zeros(length(listOfParams),2);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    
    % From tables
    A = seriesRead3ChTable(selParams,'VaryAPS',varyAPS);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    delay(iParams,:) = A.Delay;
    [delay(iParams,:),I] = sort(delay(iParams,:)); % Sort for Delays
    Yr(iParams,:) = Radii;
    Yt(iParams,:) = Thicknesses;
    discAmpl(iParams,:) = A.discAmpl;
    discAmpl(iParams,:) = discAmpl(iParams,I);
    discMeanVar(iParams,:) = A.discMeanVar;
    discMeanVar(iParams,:) = discMeanVar(iParams,I);
    discN(iParams,:) = A.discN;
    discN(iParams,:) = discN(iParams,I);
    g2vals(iParams,:) = A.g2;
    g2vals(iParams,:) = g2vals(iParams,I);
    g2std(iParams,:) = A.g2std;
    g2std(iParams,:) = g2std(iParams,I);
    
    % From DelayMeanVarX Plots
    if strcmp(typestr,'MeanVarSigma') || ...
            strcmp(typestr,'ThicknessMeanVarSigma')
        selstr = selParamsToStr(selParams);
        filelist = dir([figurepath,'*-DelayMeanVarX-',selstr,'.fig']);
        filelist = {filelist.name};
        fig = openfig([figurepath,filelist{1}]);
        figData = get(gca,'Children');
        fitStr = strjoin(figData(1).String);
        toks = regexpi(fitStr,'s =\s*([\d.-]*)','tokens');
        sigmas(iParams) = str2double(cell2mat(toks{1}));
        toks = regexpi(fitStr,'s =[\d.\s]*\(([\d.-]*)','tokens');
        sigmaConf(iParams,1) = str2double(cell2mat(toks{1}));
        toks = regexpi(fitStr,'s =[\d.\s]*\([\d.]*,\s([\d.-]*)','tokens');
        sigmaConf(iParams,2) = str2double(cell2mat(toks{1}));
        close(fig);
    end
    
    % From DelayDiscAmpl Plots
end
[~,I] = sort(delay(:,1)); % Sort for Radii
delay = delay(I,:); %delays
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
discAmpl = discAmpl(I,:);
discMeanVar = discMeanVar(I,:);
discN = discN(I,:);
g2vals = g2vals(I,:);
g2std = g2std(I,:);

%% Create figure
fig = figure;
%formatFigA5(fig);
switch typestr
    case 'Amplitude'
        for i = 1:length(I)
            plot(delay(i,:),discAmpl(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','northwest');
        xlabel(['Delay (' xUnit ')']); 
        ylabel('Coherent Amplitude'); 
        if ~varyAPS
            title('Coherent Amplitude vs. Radius of Postselected Fullcircle');
        else
            title('Coherent Amplitude vs. A_c');
        end
    case 'MeanVar'
        for i = 1:length(I)
            plot(delay(i,:),discMeanVar(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','southeast');
        ylabel('Average Variance');
        xlabel(['Delay (' xUnit ')']);
        if ~varyAPS
            title('Variance vs. Radius of Postselected Fullcircle');
        else
            title('Variance vs. A_c');
        end
    case 'MeanVarSigma'
        errorbar(Yr(:,1),sigmas,abs(sigmas-sigmaConf(:,1)), ...
            abs(sigmas-sigmaConf(:,2)),'o-','DisplayName', ...
            'Standard Deviation of Gaussian with 95% confidence intervals');
        xlabel('Ring Radius');
        ylabel('Temporal Width of Minimum Variance');
        title('Width of Minimum Variance vs. Postselected Radius');
        legend('show');
    case 'DiscN'
        for i = 1:length(I)
            plot(delay(i,:),discN(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','northwest');
        xlabel(['Delay (' xUnit ')']);
        ylabel('Photon Number');
        if ~varyAPS
            title('Photon Number vs. Radius of Postselected Fullcircle');
        else
            title('Photon Number vs. A_c');
        end
    case 'G2'
        for i = 1:length(I)
            plot(delay(i,:),g2vals(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','northeast');
        xlabel(['Delay (' xUnit ')']);
        ylabel('g^{(2)}(0)');
        if ~varyAPS
            title('G2 vs. Radius of Postselected Fullcircle');
        else
            title('G2 vs. A_c');
        end
    case 'ThicknessMeanVar'
        surf(delay,Yt,discMeanVar);
        view(-50,20);
        xlabel(['Delay (' xUnit ')']);
        zlabel('Average Variance');
        title('Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessMeanVarSigma'
        errorbar(Yt(:,1),sigmas,abs(sigmas-sigmaConf(:,1)), ...
            abs(sigmas-sigmaConf(:,2)),'o','DisplayName', ...
            'Standard Deviation of Gaussian with 95% confidence intervals');
        xlabel('Ring Thickness');
        ylabel('Temporal Width of Minimum Variance');
        title('Width of Minimum Variance vs. Postselected Thickness');
        legend('show');
    case 'ThicknessMeanVarMin'
        [~,iMin] = min(discMeanVar(1,:));
        plot(Yt(:,iMin),discMeanVar(:,iMin));
        xlabel('Ring Thickness');
        ylabel('Minimum of Postselected Variance');
        title('Minimum Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessG2'
        surf(delay,Yt,g2vals);
        view(3);
        xlabel(['Delay (' xUnit ')']);
        zlabel('g^{(2)}');
        title('Photon Number vs. Thickness of Postselected Fullcircle');
end
set(fig,'Color','w');
graphicsSettings;
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    close all;
end

end

