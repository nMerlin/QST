

chAssign = [2,1,3];

% if exist('X1','var')
%         % make sure all have the same number of pulses. 
%          if ~isequal(size(X1,1),size(X2,1),size(X3,1))
%              X1 = X1(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
%              X2 = X2(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
%              X3 = X3(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
%          end
%  
%         % set which channel ist the target channel etc
%         quadratures = zeros([size(X1) 3]);
%         quadratures(:,:,:,1) = X1;
%         quadratures(:,:,:,2) = X2;
%         quadratures(:,:,:,3) = X3;
% %         Xtg = quadratures(:,:,:,chAssign(1));  
% %         XpsFast = quadratures(:,:,:,chAssign(2));
% %         XpsSlow = quadratures(:,:,:,chAssign(3));
%         clear('quadratures'); 
% end
%%
% phaseselection = -pi:0.5:pi ;
phaseselection = 0:0.5:2*pi ;

% radius and thickness for amplitude selection
r = 10;
d = 2;
w = 0.5;
% %normalisation of the Quadratures from postselection channels
% O1norm = O1*sqrt(nPsMean/nPsFast);
% O2norm = O2*sqrt(nPsMean/nPsSlow);

%nTg = nPhotons(X2,X2,X2);
[nDs,nDsStd] = deal(length(phaseselection));
for i  = 1:length(phaseselection)  

       selParams = struct('Type','phaseAndAmplitude','Position',[phaseselection(i),w,r,d]);
       %selParams = struct('Type','phaseAndAmplitude','Position',[0.5,0.01,r,d]);
  % [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);%,'Plot','show','Filename',['test' num2str(phaseselection(i)) '-assessTheta']);
    %   [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams,'Plot','show');
     [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams);
        thetaMiraSel = oThetaMira(iSelect);
       % thetaMiraSel = selTheta;
            
     % compute photon number of postselected quadratures in the doubleslit
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
    
    nDs(i) = mean(nValues);
    nDsStd(i) = std(nValues);
end
%%
%Selection of 'fullcircle' in selectRegion for getting a mean Photonnumber
%nDs_amplitudeSelection 
selParams2 = struct('Type','fullcircle','Position',[r,d]);
      % selParams = struct('Type','phase','Position',[phaseselection(i),0.5]);
%         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams); %,'Plot','show','Filename',['test' num2str(phaseselection(i)) '-assessTheta']);
%       
        [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams2);
       thetaMiraSel = oThetaMira(iSelect);
        % thetaMiraSel = selTheta;
     % compute photon number of postselected quadratures in the doubleslit
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
    
    nDs_amplSel = mean(nValues);
    nDsStd_amplSel = std(nValues);    
    
    %%


figure(20);
% %hold on;
errorbar(phaseselection,nDs./nDs_amplSel,nDsStd,'o-','DisplayName',['-r-' num2str(r) '-d-' num2str(d) '-w-' num2str(w)]);
%errorbar(phaseselection,nDs./nDs_amplSel,nDsStd,'o-');
 %ylim([0.7 1.5]);
 xlabel('postselected phase');
 ylabel('postselected photon number');
 legend();
 grid on;
% title(['Position ',num2str(positions(j)),'mm']);
% savefig([num2str(positions(j)) 'mm-photonNumberVsPhase.fig']);
% print([num2str(positions(j)),'mm-photonNumberVsPhase.png'],'-dpng','-r300');
% clf();
        

%fit the phaseCheckPlotting
% figure(16);
% hold on;
% %fit with a sine
% sinFkt = 'a*sin(x-b)-c';
% startPoints=[0.1 0 1];
% nDsNorm=nDs./nDs_amplSel;
% f=fit(phaseselection', nDsNorm',sinFkt,'Lower',[0,0,-3],'Upper',[10,2*pi,3], 'Start', startPoints);
% plot(f,phaseselection,nDs./nDs_amplSel,'o-');


%Plots for variating the postselection range of phase and amplitude

grid on;
%title('phasewidth 0.5  amplitude r=14, d=1');
savefig('9,500mm_p0.5_ar14_ad1.fig');