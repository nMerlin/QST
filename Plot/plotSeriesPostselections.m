function plotSeriesPostselections(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultType = 'Amplitude';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[typestr] = c{:};

%% Gather data
[X,Y,Z] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    A = seriesRead3ChTable(selParams);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    X(iParams,:) = A.Delay;
    [X(iParams,:),I] = sort(X(iParams,:)); % Sort for Delays
    Y(iParams,:) = Radii;
    Z(iParams,:) = A.discAmpl;
    Z(iParams,:) = Z(iParams,I);
end
[~,I] = sort(X(:,1)); % Sort for Radii
X = X(I,:);
Y = Y(I,:);
Z = Z(I,:);

%% Create plots
switch typestr
    case 'Amplitude'
        waterfall(X,Y,Z);
        view(-20,20);
end

end

