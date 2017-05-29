function [ expQ, expP, expQ2, expP2, expDelQ, expDelP, unc, nPhotons, meanN ] = computeExpectations2( X, theta, filename )
%COMPUTEEXPECTATIONS Computes <Q>, <P>, <Q^2>, <P^2> and uncertainties,
%plots them over theta.
%
%   X and THETA should be discretized for theta with DISCRETIZETHETA
%   PHASEQ gives the arbitrary phase for Q

Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 

[~, nIntervals, nSegments] = size(X);

% Compute <Q> and <P>
meanX = zeros(nIntervals, nSegments);
meanX2 = meanX;
meanTheta = zeros(nIntervals, nSegments);
expQ = zeros(nIntervals, nSegments);
expP = expQ;
expQ2 = expQ;
expP2 = expQ;
expDelQ = expQ;
expDelP = expQ;
for iSegment = 1 : nSegments
    % Compute averages and quadratic averages
    for iInterval = 1 : nIntervals
        meanX(iInterval, iSegment) = ...
            mean(X(~isnan(X(:,iInterval,iSegment)), iInterval, iSegment));
        meanTheta(iInterval, iSegment) = ...
            mean(theta(~isnan(theta(:, iInterval, iSegment)), ...
            iInterval, iSegment));
        meanX2(iInterval, iSegment) = ...
            mean(X(~isnan(X(:,iInterval,iSegment)), ...
            iInterval, iSegment).^2);
    end % iInterval
    
    % Extract <Q>, <P>, <Q^2>, <P^2>
    expQ(:,iSegment) = meanX(:,iSegment);
    expQ2(:,iSegment) = meanX2(:,iSegment);
    expDelQ(:,iSegment) = sqrt(expQ2(:,iSegment)-(expQ(:,iSegment)).^2);  
    [~, indexP] = min(abs(meanTheta(:,iSegment) - (meanTheta(1,iSegment) + pi/2)));
    expP(:,iSegment) = [meanX(indexP:end, iSegment)' meanX(1:indexP-1, iSegment)']';
    expP2(:,iSegment) = [meanX2(indexP:end, iSegment)' meanX2(1:indexP-1, iSegment)']';
    expDelP(:,iSegment) = sqrt(expP2(:,iSegment)-(expP(:,iSegment)).^2);
end % iSegment 

% averages meanX over segments, then extracts expectations
% for iSegment = 1 : nSegments
%     % Compute averages and quadratic averages
%     for iInterval = 1 : nIntervals
%         meanX(iInterval, iSegment) = ...
%             mean(X(~isnan(X(:,iInterval,iSegment)), iInterval, iSegment));
%         meanTheta(iInterval, iSegment) = ...
%             mean(theta(~isnan(theta(:, iInterval, iSegment)), ...
%             iInterval, iSegment));
%         meanX2(iInterval, iSegment) = ...
%             mean(X(~isnan(X(:,iInterval,iSegment)), ...
%             iInterval, iSegment).^2);
%     end % iInterval
% end % iSegment
% meanX=mean(meanX,2);
% meanX2=mean(meanX2,2);
% meanTheta=mean(meanTheta,2);
% expQ = meanX;  % Extract <Q>, <P>, <Q^2>, <P^2>
% expQ2 = meanX2;
% expDelQ = sqrt(expQ2-(expQ).^2);  
% [~, indexP] = min(abs(meanTheta - (meanTheta(1) + pi/2)));
% expP = [meanX(indexP:end)' meanX(1:indexP-1)']';
% expP2 = [meanX2(indexP:end)' meanX2(1:indexP-1)']';
% expDelP = sqrt(expP2-(expP).^2);

%Segment to be used
Seg = 1;

% Coherent state photon number
nPhotons = (expQ2 + expP2)/(4*Norm^2) -0.5 ;
meanN = mean(nPhotons(:,Seg));
%meanN = mean(nPhotons);  %if averaged over segments

%uncertainty value
unc = expDelQ.*expDelP;

% Plot
close all;
x = meanTheta(:,Seg);
plot(x, expP(:,Seg), x, expP2(:,Seg), x, expQ(:,Seg), x, expQ2(:,Seg), x, expDelQ(:,Seg),...
x, expDelP(:,Seg), x,  unc(:,Seg), x, nPhotons(:,Seg), 'linewidth', 2);
% x = meanTheta;
% plot(x, expP, x, expP2, x, expQ, x, expQ2, x, expDelQ,...
% x, expDelP, x,  unc, x, nPhotons, 'linewidth', 2);   %if averaged over segments

hold on;
plot(x,Norm^2*ones(length(x)),'k-','lineWidth',0.5);
xlabel('Phase \theta');
axis([0 2*pi min(min(meanX))-0.5 max(max(meanX2))+0.5]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<P>$', '$<P^{2}>$', '$<Q>$', '$<Q^{2}>$', '$\Delta Q$',...
    '$\Delta P$', '$\Delta Q \cdot \Delta P$ ', '$<n>$ ','uncertainty limit', 'location', 'best');
title(strcat('Mean Photon number =','\, ',num2str(meanN)));

% Save plot
print(strcat('expect-',filename,'-', strrep(num2str(meanN),'.',','),'photons'), '-dpng');

end % computeExpectations

