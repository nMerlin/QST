function [ expQ, expP, expQ2, expP2, expDelQ, expDelP, unc, nPhotons, meanN ] = computeExpectations2( X, theta )
%COMPUTEEXPECTATIONS Computes <Q>, <P>, <Q^2>, <P^2>
%
%   X and THETA should be discretized for theta with DISCRETIZETHETA
%   PHASEQ gives the arbitrary phase for Q

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

% Coherent state photon number
nPhotons = expQ.^2 + expP.^2;
meanN = mean(mean(nPhotons));

%uncertainty value
unc = expDelQ.*expDelP;

% Plot
for i = 1:1 %nSegments
     x = meanTheta(:,i);
     plot(x, expP(:,i), x, expP2(:,i), x, expQ(:,i), x, expQ2(:,i), x, expDelQ(:,i),...
     x, expDelP(:,i), x,  unc(:,i), x, nPhotons(:,i), 'linewidth', 2);
     xlabel('Phase \theta');
end
axis([0 2*pi min(min(meanX))-0.5 max(max(meanX2))+0.5]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<P>$', '$<P^{2}>$', '$<Q>$', '$<Q^{2}>$', '$\Delta Q$',...
    '$\Delta P$', '$\Delta Q \cdot \Delta P$ ', '$<n>$ ', 'location', 'best');
title(strcat('Mean Photon number =','\, ',num2str(meanN)));

% Save plot
%dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
%print(['expectations-' dateString], '-dpng');

end % computeExpectations

