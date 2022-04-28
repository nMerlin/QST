function[thetaCorr,thetaMiraCorr,XtgCorr,XpsSlowCorr,XpsFastCorr]=corrBoundariesOfPiezosegments(Xtg,XpsSlow,XpsFast,piezoSign);
%function[XtgCorr,XpsSlowCorr,XpsFastCorr]=corrBoundariesOfPiezosegments(Xtg,XpsSlow,XpsFast,piezoSign);

%the theta's and X are flipped in timing sequenz in every second
%PiezoSegment

[a,b,c]=size(Xtg);
seg = a*b;
piezoSeg=c;
% first = 102952;
% last = 20800;
first = 102952;
last = 20800;
[XtgNew] = flipPiezoSegments(Xtg);
[XpsSlowNew] = flipPiezoSegments(XpsSlow);
[XpsFastNew] = flipPiezoSegments(XpsFast);


%%
%theta and thetaMira are calculated by using computePhase with the flipped
%X's
periodsPerSeg = 2;

[thetaNew,~] = computePhase(XtgNew,XpsSlowNew,piezoSign,'Period',periodsPerSeg); 

[thetaMiraNew,~] = computePhase(XtgNew,1,piezoSign,'Period',periodsPerSeg,'Peakthreshold',0.1);

%thetaNew and thetaMiraNew are corrected by cutting the boundaries for each
%PiezoSegment
thetaCorr=thetaNew;
thetaCorr(first:seg,:)=[];
thetaCorr(1:last,:)=[];

thetaMiraCorr=thetaMiraNew;
thetaMiraCorr(first:seg,:)=[];
thetaMiraCorr(1:last,:)=[];
%%
%the Quadratures X are corrected by cutting the boundaries for each
%PiezoSegment
XtgNew = reshape(XtgNew,[seg,piezoSeg]);
XtgCorr=XtgNew;
XtgCorr(first:seg,:)=[];
XtgCorr(1:last,:)=[];

XpsSlowNew = reshape(XpsSlowNew,[seg,piezoSeg]);
XpsSlowCorr=XpsSlowNew;
XpsSlowCorr(first:seg,:)=[];
XpsSlowCorr(1:last,:)=[];

XpsFastNew = reshape(XpsFastNew,[seg,piezoSeg]);
XpsFastCorr=XpsFastNew;
XpsFastCorr(first:seg,:)=[];
XpsFastCorr(1:last,:)=[];
%%
% %find if thetaCorr is more often bigger than thetaMiraCorr or vice versa
% %for making plots with theta-thetaMira better - for the calculation of
% %photonnumber with uniformSampling.m this correction causes problems
% md =  median(thetaCorr-thetaMiraCorr);
% %correct for cases where thetaCorr and thetaMiraCorr are falsely in the
% %wrong order due to mod2pi.
% for ps =1:piezoSeg
%     if md(ps) > 0
%         thetaCorr(thetaMiraCorr(:,ps)>thetaCorr(:,ps),ps)=thetaCorr(thetaMiraCorr(:,ps)>thetaCorr(:,ps),ps)+2*pi;
%     else
%         thetaMiraCorr(thetaCorr(:,ps)>thetaMiraCorr(:,ps),ps)=thetaMiraCorr(thetaCorr(:,ps)>thetaMiraCorr(:,ps),ps)+2*pi;
%     end
% end
end