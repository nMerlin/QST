function [ expQ, expP, expQ2, expP2, expDelQ, expDelP, unc, nPhotons ] = computeExpectations( X, theta, phaseQ )
%COMPUTEEXPECTATIONS Computes <Q>, <P>, <Q^2>, <P^2>
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
expQ = zeros(1, nSegments);
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
    [~, indexQ] = min(abs(meanTheta(:,iSegment) - phaseQ));
    expQ(iSegment) = meanX(indexQ, iSegment);
    expQ2(iSegment) = meanX2(indexQ, iSegment);
    expDelQ(iSegment) = sqrt(expQ2(iSegment)-(expQ(iSegment))^2);
    [~, indexP] = min(abs(meanTheta(:,iSegment) - (phaseQ + pi/2)));
    expP(iSegment) = meanX(indexP, iSegment);
    expP2(iSegment) = meanX2(indexP, iSegment);
    expDelP(iSegment) = sqrt(expP2(iSegment)-(expP(iSegment))^2);
end % iSegment

% Coherent state photon number
nPhotons = (expQ2 + expP2)/(4*Norm^2) -0.5 ;

%uncertainty value
unc = expDelQ.*expDelP;

% Plot
x = 1 : nSegments;
plot(x, expP, x, expP2, x, expQ, x, expQ2, x, expDelQ, x, expDelP, x,  unc, x, nPhotons, 'linewidth', 2);
xlabel('Piezo Segments');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<P>$', '$<P^{2}>$', '$<Q>$', '$<Q^{2}>$', '$\Delta Q$',...
    '$\Delta P$', '$\Delta Q \cdot \Delta P$ ', 'photons', 'location', 'best');

% Save plot
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
print(['expectations-' dateString], '-dpng');

end % computeExpectations

