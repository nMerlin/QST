t = 0.4;
listOfParams = struct('Type',{'fullcircle'},'Position',{[0.25 t],[0.5 t],[1,t],[1.5 t],...
[2 t],[2.5 t],[3 t],[3.5 t],[4,t],[5 t],[6 t],[7 t],[8 t],[9 t],[10 t],[11 t],[12 t],[13 t]});
discAmpl = zeros(length(listOfParams),1);
[Aps,varQ,varP,discN,meang2] = deal(zeros(length(listOfParams),1));
nDisc = 100;
chAssign = [3,1,2]; %[target,ps_piezo_fast,ps_piezo_slow]
 if exist('X1','var')
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1)); 
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
 end

 if  (~exist('theta','var'))
     [theta,~] = computePhase(Xtg,XpsFast,piezoSign);
 end
 if (~exist('O1','var'))   %%case, where remMod is done first?
            [O1,O2,O3,oTheta,iOrth] = selectOrthogonal(XpsFast,XpsSlow,Xtg,theta,piezoSign);
 end
 
if remMod
    n1vec = photonNumberVector(X1);
    n2vec = photonNumberVector(X2);
    n3vec = photonNumberVector(X3);
    n1vec = n1vec(iOrth);
    n2vec = n2vec(iOrth);
    n3vec = n3vec(iOrth);
    [O3rem,O1rem,O2rem,oThetaRem,iSelect] = removeNBelowLimit(O3,O1,O2,oTheta,n1vec,n2vec,n3vec,0.5);
end
for i = 1:length(listOfParams)
    if remMod
        [selX,selTheta] = selectRegion(O1rem,O2rem,O3rem,oThetaRem,listOfParams(i));
    else
        [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,listOfParams(i));
    end
    %[selX,selTheta] = selectRegion(O1remMod,O2remMod,O3remMod,oThetaRemMod,listOfParams(i));
    [disSelX,disSelTheta]=discretizeTheta(selX,selTheta,nDisc);
        [expQ,expP,expQ2,expP2,delQ,delP,~,~,~,~] = ...
            computeExpectations2(disSelX,disSelTheta,'bla','Plot','hide');
    discAmpl(i) = expQ(round(nDisc/4));
    selParams = listOfParams(i);
    Aps(i) = selParams.Position(1);
    varVector =  expQ2(:)-(expQ(:)).^2; 
    varQ(i) = mean(varVector(round(nDisc/4)-2:round(nDisc/4)+2));
    varP(i) = mean(varVector(round(nDisc/2)-2:round(nDisc/2)+2));
    % g2 estimation
    [g2values,nValues] = deal(zeros(10,1));
    for iG2=1:10
        try
            uniformX = seriesUniformSampling(selX,selTheta,'NBins',100);
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling for g2.']);
            uniformX = selX;
        end    
        try
            [g2values(iG2),nValues(iG2)] = g2(uniformX,length(uniformX));
        catch
            warning('not enough X to compute g2. Set g2 = 2.');
            g2values(iG2) = 2;
            uniformX=uniformX-mean(uniformX);
            nValues(iG2) = mean(uniformX.^2)-0.5;
        end
    end
    
    discN(i) = mean(nValues);
    meang2(i) = mean(g2values);
end

plot(Aps,discAmpl,'o-');
xlabel('A_{ps}');ylabel('coherent Amplitude');
graphicsSettings;
savefig([filename '-thickness-' num2str(t) '-Ampl.fig']);
print([filename '-thickness-' num2str(t) '-Ampl.png'],'-dpng');
clf();

plot(Aps,discN,'o-');
xlabel('A_{ps}');ylabel('postselected Photon Number');
graphicsSettings;
savefig([filename '-thickness-' num2str(t) '-n.fig']);
print([filename '-n.png'],'-dpng');
clf();

plot(Aps,varQ,'o-',Aps,varP,'o-');
xlabel('A_{ps}');ylabel('Variance');
legend('var_{Q}','var_{P}');
graphicsSettings;
savefig([filename '-thickness-' num2str(t) '-Var.fig']);
print([filename '-thickness-' num2str(t) '-Var.png'],'-dpng');
clf();

plot(Aps,meang2,'o-');
xlabel('A_{ps}');ylabel('g^{(2)}(0)');
graphicsSettings;
savefig([filename '-thickness-' num2str(t) '-g2.fig']);
print([filename '-thickness-' num2str(t) '-g2.png'],'-dpng');
clf();


