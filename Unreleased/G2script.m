clear all;
[X1, X2, ~, piezoSign, ~] = prepareData('01-56,985mm-5mW-LOonly.raw','02-56,985mm-5mW-LOwithDL.raw','Channels',1:2,'Offset','global','Piezo','yes');

offset1 = (mean(max(max(max(X1,1))))+mean(min(min(min(X1,1)))))/2;
X1 = X1 - offset1;
offset2 = (mean(max(max(max(X2,1))))+mean(min(min(min(X2,1)))))/2;
X2 = X2 - offset2;

[theta1,~,~] = computePhase(X1,ones(size(X1)),piezoSign,'Period',1);
[theta2,~,~] = computePhase(X2,ones(size(X2)),piezoSign,'Period',2);
% [theta1,theta10,~] = computePhase(X1,ones(size(X1)),piezoSign,'Period',2);
% [theta2,theta20,~] = computePhase(X2,ones(size(X2)),piezoSign,'Period',4);
[uniformX1,uniformTheta1] = seriesUniformSampling(X1(:),theta1(:),'NBins',200);
[uniformX2,uniformTheta2] = seriesUniformSampling(X2(:),theta2(:),'NBins',200);
theta12 = mod(theta2-theta1,2*pi);
[uniformX12,uniformTheta12] = seriesUniformSampling(X1(:)+X2(:),theta12(:),'NBins',300);


[g2values,nPhotons] = g2(uniformX1,length(uniformX1))
[g2values2,nPhotons2] = g2(uniformX2,length(uniformX2))
g2values12 = g22Ch(uniformX1, uniformX2, uniformX12,length(uniformX1),length(uniformX2),length(uniformX12))

% theta120 = mod(theta20-theta10,2*pi);
% [uniformX120,uniformTheta120] = seriesUniformSampling(X1(:)+X2(:),theta120(:),'NBins',300);
% g2values120 = g22Ch(uniformX1, uniformX2, uniformX120,length(uniformX1),length(uniformX2),length(uniformX120));
% 
% histogram(uniformTheta12);
% plot(uniformTheta12,uniformX12(:));
% histogram(uniformTheta1);  
% histogram(uniformTheta2);

% X1Start = X1(:);
% X2Start = X2(:);
% theta1Start = theta1(:);
% theta2Start = theta2(:);
% %theta12Start = theta12(:);
% theta12Start =  mod(theta2(:)-theta1(:),2*pi);
% iterN = 3;
% for i = 1: iterN
% [~,~,uniformIndicesMidX12] = seriesUniformSamplingIndex(X1Start+X2Start,theta12Start,'NBins',300);
% X1Mid = X1Start(uniformIndicesMidX12);
% X2Mid = X2Start(uniformIndicesMidX12);
% theta1Mid = theta1Start(uniformIndicesMidX12);
% theta2Mid = theta2Start(uniformIndicesMidX12);
% 
% [X1End,theta1End, indicesX1End] = seriesUniformSamplingIndex(X1Mid,theta1Mid,'NBins',300);
% X1Mid = X1End;
% theta1Mid = theta1End;
% X2Start = X2Mid(indicesX1End);
% theta2Start = theta2Mid(indicesX1End);
% 
% [X2End,theta2End, indicesX2End] = seriesUniformSamplingIndex(X2Start,theta2Start,'NBins',300);
% X1End = X1Mid(indicesX2End);
% theta1End = theta1Mid(indicesX2End);
% 
% theta12Start = mod(theta2End-theta1End,2*pi);
% X1Start = X1End;
% X2Start = X2End;
% theta1Start = theta1End(:);
% theta2Start = theta2End(:);
% %theta12Start = theta12End(:);
% %[uniformX12subsub,uniformTheta12subsub,~] = seriesUniformSamplingIndex(X1End+X2subsub,theta12subsub,'NBins',300);
% end
tic;
[uniform2ChX1,uniform2ChX2,uniform2ChTheta1,uniform2ChTheta2] = uniformSampling2Channel(X1,X2,theta1,theta2,5);
theta2Ch12 =  mod(uniform2ChTheta1-uniform2ChTheta2,2*pi);
[uniform2ChX12,uniform2ChTheta12,~] = seriesUniformSamplingIndex(uniform2ChX1+uniform2ChX2,theta2Ch12,'NBins',100);
g2values12sum = g22Ch(uniform2ChX1, uniform2ChX2, uniform2ChX1+uniform2ChX2,length(uniform2ChX1),length(uniform2ChX2),length(uniform2ChX1))
g2values12unif = g22Ch(uniform2ChX1, uniform2ChX2,uniform2ChX12,length(uniform2ChX1),length(uniform2ChX2),length(uniform2ChX12))
[g2values2Ch1,nPhotons2Ch1] = g2(uniform2ChX1,length(uniform2ChX1))
[g2values2Ch2,nPhotons2Ch2] = g2(uniform2ChX2,length(uniform2ChX2))
toc;
%[uniformX12,uniformTheta12,~] = seriesUniformSamplingIndex(X1Start+X2Start,theta12Start,'NBins',300);

%[uniformX2subsub,uniformTheta2subsub,~] = seriesUniformSamplingIndex(X2subsub,theta2subsub,'NBins',300);

%g2values12sum = g22Ch(X1Start, X2Start, X1Start+X2Start,length(X1Start),length(X2Start),length(X1Start+X2Start));
%g2values12unif = g22Ch(X1Start, X2Start,uniformX12,length(X1Start),length(X2Start),length(uniformX12));



