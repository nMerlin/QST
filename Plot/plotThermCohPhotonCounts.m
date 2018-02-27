function [thermal,coherent] = plotThermCohPhotonCounts(n,nAv)
%PLOTTHERMCOHPHOTONCOUNTS Creates bar plots of thermal and coherent photon
%counts.
%
%   Input Arguments:
%       n - Number samples to simulate
%       nAv - Average photon count
%
%   Output Arguments:
%       thermal - vector containing thermal photon counts
%       coherent - vector containing coherent photon counts

threshold = 0.001;
thermal = zeros(n,1);
coherent = zeros(n,1);
g2 = 0;

%% Drawing customized samples
while abs((mean(thermal)-nAv)/nAv)>threshold || abs(g2-2)>0.01
    thermal = exponentialCDF(rand(n,1),nAv,'Inverse',true);
    g2 = mean(thermal.^2)/mean(thermal)^2;
end
while abs((mean(coherent)-nAv)/nAv)>threshold
    coherent = normrnd(nAv,sqrt(nAv),n,1);
end

%% Creating plots
close all
figTh = figure;
figTh.Units = 'pixels'; 
figTh.Position = [230,600,800,200];
mTh = mean(thermal);
bTh = bar(thermal,'y');
bTh.FaceColor = [246/255,189/255,22/255];
set(gca,'XLim',[0 n+1]);
set(gca,'YLim',[0 max(thermal)]);
set(gca,'Box','off','XTick',[],'YColor','None','LineWidth',2);
line([1 n],[mTh mTh],'LineStyle',':','Color','black','LineWidth',2);
legend('Photon Counts','Average');
figCo = figure;
figCo.Units = 'pixels'; 
figCo.Position = [230,250,800,200];
mCo = mean(coherent);
bCo = bar(coherent);
bCo.FaceColor = [132/255,184/255,25/255];
set(gca,'XLim',[0 n+1]);
set(gca,'YLim',[0 max(thermal)]);
set(gca,'Box','off','XTick',[],'YColor','None','LineWidth',2);
line([1 n],[mCo mCo],'LineStyle',':','Color','black','LineWidth',2);
legend('Photon Counts','Average');

end

