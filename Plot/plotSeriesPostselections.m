function plotSeriesPostselections(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultType = 'Amplitude';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,typestr] = c{:};

%% Gather data
[X,Yr,Yt,discAmpl,discMeanVar,discN,g2vals,g2std] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    A = seriesRead3ChTable(selParams);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    X(iParams,:) = A.Delay;
    [X(iParams,:),I] = sort(X(iParams,:)); % Sort for Delays
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
end
[~,I] = sort(X(:,1)); % Sort for Radii
X = X(I,:);
Yr = Yr(I,:);
Yt = Yt(I,:);
discAmpl = discAmpl(I,:);
discMeanVar = discMeanVar(I,:);
discN = discN(I,:);
g2vals = g2vals(I,:);
g2std = g2std(I,:);

%% Create figure
fig = figure;
formatFigA5(fig);
switch typestr
    case 'Amplitude'
        waterfall(X,Yr,discAmpl);
        view(-20,20);
        xlabel('Delay (fs)');
        zlabel('Coherent Amplitude');
        title('Coherent Amplitude vs. Radius of Postselected Fullcircle');
    case 'MeanVar'
        surf(X,Yr,discMeanVar);
        view(-50,20);
        xlabel('Delay (fs)');
        zlabel('Average Variance');
        title('Variance vs. Radius of Postselected Fullcircle');
    case 'DiscN'
        waterfall(X,Yr,discN);
        view(-20,20);
        xlabel('Delay (fs)');
        zlabel('Photon Number');
        title('Photon Number vs. Radius of Postselected Fullcircle');
    case 'G2'
        waterfall(X,Yr,g2vals);
        view(3);
        xlabel('Delay (fs)');
        zlabel('g^{(2)}');
        title('Photon Number vs. Radius of Postselected Fullcircle');
    case 'ThicknessMeanVar'
        surf(X,Yt,discMeanVar);
        view(-50,20);
        xlabel('Delay (fs)');
        zlabel('Average Variance');
        title('Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessMeanVarMin'
        [~,iMin] = min(discMeanVar(1,:));
        plot(Yt(:,iMin),discMeanVar(:,iMin));
        xlabel('Ring Thickness');
        ylabel('Minimum of Postselected Variance');
        title('Minimum Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessG2'
        waterfall(X,Yt,g2vals);
        view(3);
        xlabel('Delay (fs)');
        zlabel('g^{(2)}');
        title('Photon Number vs. Thickness of Postselected Fullcircle');
end
set(fig,'Color','w');

%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    close all;
end

end

