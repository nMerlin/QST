function [] = Expectations3ChNumberSeries()


%This script writes the expectation values for a photon number series of 
%3 Channels into an excelsheet. The data X1,X2,X3, theta and piezoSign
%should be stored in files like '2017-09-20-validated-04-5mW-LOwithDL.mat'.

date = '2017-09-20';

%% iteration over data files
number = 56:2:90; %number of data file
[nX1, nX2, nX3, meanVarX, q1, q21, p1, p21, qmax,q2max,varQ,varP] = deal(zeros(length(number),1));

for i = 1:length(number)
    if number(i) <10
        numberName = ['0' num2str(number(i))];
    else
        numberName = num2str(number(i));
    end
    
    % load data
    load([date '-validated-' numberName '-5mW-LOwithDL.mat']);
    
    % compute photon numbers for each channel
    [n1,n2,n3] = nPhotons(X1,X2,X3);
    nX1(i) = n1;
    nX2(i) = n2;
    nX3(i) = n3;
    
    % postselections
    [O1,O2,O3,oTheta] = selectOrthogonal(X1,X2,X3,theta,piezoSign);
    filename = [date '-' numberName];
    [selX,selTheta, ~, mVarX] = selectRegion(O1,O2,O3,oTheta,'Type','fullcircle',...
        'Position',[2.5 0.5],'Plot','show','Output','print','Filename',filename);
    close all;
    meanVarX(i) = mVarX;
    
    % descretization and expectation values
    [selX,selTheta]=discretizeTheta(selX,selTheta,100);
    [expQ,expP,expQ2,expP2,delQ,delP,~,~,~,~] = ...
        computeExpectations2(selX,selTheta,'bla','Plot','hide');
    q1(i) = expQ(1);
    q21(i) = expQ2(1);
    p1(i) = expP(1);
    p21(i) = expP2(1);
    qmax(i) = max(expQ);
    q2max(i) = max(expQ2);
    varQ(i) = delQ^2;
    varP(i) = delP^2;
    
end

%% write results in excel sheet
M = [number' nX1 nX2 nX3 meanVarX q1 q21 p1 p21 qmax q2max varQ varP];
xlswrite('test',M);

end