  
 theta = 0 + 2*pi*rand(10000,1);  
 intens = normrnd(0,20,[10000,1]);
 x=zeros(10000,1);
 y=zeros(10000,1);

for i=1:length(x)
    x(i) = intens(i)*cos(theta(i));
    y(i) = intens(i)*sin(theta(i));
end
 
%amplitude and phase selection
minTheta = 0;
maxTheta = 0.5;
minInt = 10;
maxInt = 14;

hustheta =atan2(x,y);
thetaSel = find(hustheta < maxTheta & hustheta > minTheta &  sqrt(x.^2+y.^2)>minInt & sqrt(x.^2+y.^2)<maxInt);

% xsel=x(thetaSel);
% ysel=y(thetaSel);
intensSel=intens(thetaSel);
thetaSel=theta(thetaSel);

        plot(x,y,'.');
        hold on;
        plot(xsel,ysel,'.');
 
allPhases = 0:0.01:2*pi;
Int = mean(Int,1);




matrix=ones(length(position),length(Sel),length(allPhases));
matrix(1,1,:)=allPhases;
matrix(1,:,1)=thetaSel;
matrix(:,1,1)=position;

matrix=ones(length(position),length(Sel),length(allPhases));
matrix(1,1,:)=allPhases;
matrix(1,:,1)=thetaSel;
matrix(:,1,1)=position;
