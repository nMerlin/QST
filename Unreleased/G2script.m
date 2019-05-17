[X1, X2, ~, piezoSign, ~] = prepareData('05-5mW-LOonly.raw','06-5mW-LOasSIG.raw','Channels',1:2,'Offset','global','Piezo','yes');
% [theta1,theta10,~] = computePhase(X1,ones(size(X1)),piezoSign,'Period',1);
% [theta2,theta20,~] = computePhase(X2,ones(size(X2)),piezoSign,'Period',2);
[theta1,theta10,~] = computePhase(X1,ones(size(X1)),piezoSign,'Period',2);
[theta2,theta20,~] = computePhase(X2,ones(size(X2)),piezoSign,'Period',4);
[uniformX1,uniformTheta1] = seriesUniformSampling(X1(:),theta1(:),'NBins',200);
[uniformX2,uniformTheta2] = seriesUniformSampling(X2(:),theta2(:),'NBins',200);
theta12 = mod(theta2-theta1,2*pi);
[uniformX12,uniformTheta12] = seriesUniformSampling(X1(:)+X2(:),theta12(:),'NBins',300);
[g2values,nPhotons] = g2(uniformX1,length(uniformX1));
[g2values2,nPhotons2] = g2(uniformX2,length(uniformX2));
g2values12 = g22Ch(uniformX1, uniformX2, uniformX12,length(uniformX1),length(uniformX2),length(uniformX12));

theta120 = mod(theta20-theta10,2*pi);
[uniformX120,uniformTheta120] = seriesUniformSampling(X1(:)+X2(:),theta120(:),'NBins',300);
g2values120 = g22Ch(uniformX1, uniformX2, uniformX120,length(uniformX1),length(uniformX2),length(uniformX120));

histogram(uniformTheta12);
plot(uniformTheta12,uniformX12(:));
histogram(uniformTheta1);
histogram(uniformTheta2);