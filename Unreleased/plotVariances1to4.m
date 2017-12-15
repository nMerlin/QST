% Improvised script to plot variances of data sets X1,..,X4 and
% theta1,...,theta3.
Norm = 1/sqrt(2);
[~,~,~,T1]=assessTheta(theta1(:,1:3),X1(:,:,1:3),'Output','none');
[~,~,~,T2]=assessTheta(theta2(:,1:3),X2(:,:,1:3),'Output','none');
[~,~,~,T3]=assessTheta(theta3(:,1:3),X3(:,:,1:3),'Output','none');
[~,~,~,T4]=assessTheta(theta4(:,1:3),X4(:,:,1:3),'Output','none');
[~,~,~,T5]=assessTheta(theta5(:,1:3),X5(:,:,1:3),'Output','none');
xAxis = table2array(T1(:,'Phase'));
varX1 = table2array(T1(:,'varX'));
varX2 = table2array(T2(:,'varX'));
varX3 = table2array(T3(:,'varX'));
varX4 = table2array(T4(:,'varX'));
varX5 = table2array(T5(:,'varX'));
hold on;
plot(xAxis,varX1,'-o');
plot(xAxis,varX2,'-o');
plot(xAxis,varX3,'-o');
plot(xAxis,varX4,'-o');
plot(xAxis,varX5,'-o');
plot(xAxis,Norm^2*ones(size(xAxis)),'k.','lineWidth',0.5);
hold off;
legend('Var(X) of Channel 1 (17/09/08)', ...
    'Var(X) of Channel 2 (17/09/08)', ...
    'Var(X) of Channel 3 (17/09/08)', ...
    'Var(X) of Channel 1 (17/05/29)', ...
    ['Var(Coherent State ) = ',num2str(Norm^2)]);
xlabel('\theta');
title(['Phase-Binned Quadrature Values (',num2str(length(xAxis)),' Bins)']);
set(gca,'XLim',[min(xAxis) max(xAxis)]);