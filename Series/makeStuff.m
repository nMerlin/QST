tic
%% create new mat data and 3 Ch tables
%seriesQuadratures('Channels',1:3,'Offset','local','Piezo','yes');
%g2SeriesFromQuadratures(1000,'Weight','yes','UseX','X1','Parameter','delay');
%listOfParams = struct('Type',{'fullcircle'},'Position',{[1 0.25]});
% makeSelectionPlots('table','ListOfParams',listOfParams,'RecomputeTheta',false,...
%         'RecomputeOrth',true,'saveOrth',true,...
%         'SavePostselection',true,'SaveTheta',true,'GetDelay',true,'RemoveModulation',false,'XUnit','ps','VaryAPS',true);


%% Radius Series
%   listOfParams = struct('Type',{'fullcircle'},'Position',{[0.5 0.25],[1 0.25],[1,0.5],...
%         [2 0.5],[3 0.5],[4 0.5],[5 0.5],[6 0.5],[7 0.5],[8 0.5]});
%    makeSelectionPlots('radiusPlots','ListOfParams',listOfParams,'RecomputeTheta',false,...
%         'SavePostselection',false,'SaveTheta',false,'GetDelay',true,'RemoveModulation',false,'XUnit','ps','VaryAPS',true);
    
%% Thickness Series
listOfParamsT = struct('Type',{'fullcircle'},'Position',{[0 0.05],[0 0.1],[0,0.15],...
        [0 0.2],[0 0.4],[0 0.6],[0 0.8],[0 1]});    
makeSelectionPlots('radiusPlots','ListOfParams',listOfParamsT,'RecomputeTheta',false,...
    'RecomputeOrth',false,'saveOrth',false,...
    'SavePostselection',false,'SaveTheta',false,'GetDelay',true,'RemoveModulation',false,'XUnit','ps','VaryAPS',true);

toc