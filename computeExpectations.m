function [ expQ, expP, expQ2, expP2 ] = computeExpectations( X, theta, phaseQ )
%COMPUTEEXPECTATIONS Computes <Q>, <P>, <Q^2>, <P^2>
%
%   X and THETA should be discretized for theta with DISCRETIZETHETA
%   PHASEQ gives the arbitrary phase for Q

[~, nIntervals, nSegments] = size(X);

% Compute <Q> and <P>
meanX = zeros(nIntervals, nSegments);
meanX2 = meanX;
meanTheta = zeros(nIntervals, nSegments);
expQ = zeros(1, nSegments);
expP = expQ;
expQ2 = expQ;
expP2 = expQ;
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
    [~, indexQ] = min(abs(meanTheta(:,iSegment) - phaseQ));
    expQ(iSegment) = meanX(indexQ, iSegment);
    expQ2(iSegment) = meanX2(indexQ, iSegment);
    [~, indexP] = min(abs(meanTheta(:,iSegment) - (phaseQ + pi/2)));
    expP(iSegment) = meanX(indexP, iSegment);
    expP2(iSegment) = meanX2(indexP, iSegment);
end % iSegment

% Plot
x = 1 : nSegments;
plot(x, expP, x, expP2, x, expQ, x, expQ2, 'linewidth', 2);
xlabel('Piezo Segments');
legend('<P>', '<P^2>', '<Q>', '<Q^2>', 'location', 'best');

% Save plot
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
print(['expectations-' dateString], '-dpng');

end % computeExpectations

