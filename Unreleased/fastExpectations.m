rawDataContents = dir('raw-data');
names = {rawDataContents.name};
s = names{end};
nDL = str2num(s(1:2));
nDL = 2;
nLO = nDL-1;
[X1,X2,X3,piezoSign]=prepare3ChData(['0',num2str(nLO),'-5mW-LOonly.raw'],['0',num2str(nDL),'-5mW-LOwithDL.raw']);
[theta,ys] = computePhase(X1,X3,piezoSign);
% [O1,O2,O3,oTheta] = selectOrthogonal(X1,X2,X3,theta,piezoSign);
% [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,'Type','fullcircle','Position',[2.5 0.5],'Plot','show');
[nX1,nX2,nX3] = nPhotons(X1,X2,X3)
% ys = smoothCrossCorr(X1,X2);
% amplX1X2 = mean(max(ys)-min(ys))
% ys = smoothCrossCorr(X1,X3);
% amplX1X3 = mean(max(ys)-min(ys))
% ys = smoothCrossCorr(X2,X3);
% amplX2X3 = mean(max(ys)-min(ys))
iSelect = (X1<=2.5 & X1>=2.0);
selX = X3(iSelect);
selTheta = theta(iSelect);
assessTheta(selTheta,selX);