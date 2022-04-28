
folder='mat-data';

[filenames,~,positions]= getParametersFromFilenames('Folder',folder,'Parameter','position');
%fitCoeff = struct; 


for j = 1:length(filenames)

filename = cell2mat(filenames(j));
load([folder '\' filename], 'X1','X2','X3','O1','O2','O3','oTheta','oThetaMira');

chAssign = [2,1,3];
% 
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
r=10;
d=2;

nTg = nPhotons(X2,X2,X2);
[nDs,nDsStd] = deal(length(phaseselection));
for i  = 1:length(phaseselection)  

       selParams = struct('Type','phaseAndAmplitude','Position',[phaseselection(i),0.5,r,d]);
      % selParams = struct('Type','phase','Position',[phaseselection(i),0.5]);
%         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams); %,'Plot','show','Filename',['test' num2str(phaseselection(i)) '-assessTheta']);
%       
        [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams);
        thetaMiraSel = oThetaMira(iSelect);   
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
    %nDsStd_amplSel = std(nValues);
 
    
    %%


figure(8);
% hold on;
%fit with a sine
sinFkt = 'a^2*sin(x-b)-c';
startPoints=[0.2 0 1];

nDsNorm=nDs./nDs_amplSel;
f=fit(phaseselection', nDsNorm',sinFkt, 'Start', startPoints);
plot(f,phaseselection,nDs./nDs_amplSel,'o-');

%save fit parameters
fitCoeff.position(j) = positions(j);
fitCoeff.amplitude(j) = f.a^2;
fitCoeff.phaseShift(j) = f.b;
fitCoeff.offset(j) = f.c;

%Output in txt format
nDs_norm=nDs./nDs_amplSel;   
T=table(phaseselection', nDs_norm', nDsStd',nDs' , 'Variablenames',{'phaseselection', 'nDs_norm', 'nDsStd', 'nDs'});
writetable(T,[num2str(positions(j)) 'mm-photonNumberVsPhase.txt']);

%%%plot
% 
% errorbar(phaseselection,nDs./nDs_amplSel,nDsStd,'o-','DisplayName',['pos-' num2str(positions(j)) '-r-' num2str(r) '-d-' num2str(d)]);
% ylim([0.6 1.6]);
% xlabel('postselected phase');
% ylabel('postselected photon number');
% grid on;
% hold on;
% legend();
% title(['Position ',num2str(positions(j)),'mm']);
% savefig([num2str(positions(j)) 'mm-photonNumberVsPhase.fig']);
% print([num2str(positions(j)),'mm-photonNumberVsPhase.png'],'-dpng','-r300');
% hold off;
% clf();
end
            
