function [ expQ, expP, expQ2, expP2, delQ, delP, meanUnc, nPhotons, meanN, nSegments ] = computeExpectations2( X, theta, filename )
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
nPhotons = expQ;
unc = expQ;
meanN = zeros(nSegments);
delQ = meanN;
delP = meanN;
meanUnc = meanN;

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

    % Coherent state photon number
    nPhotons(:,iSegment) = (expQ2(:,iSegment) + expP2(:,iSegment))/(4*Norm^2) -0.5 ;
    meanN(iSegment) = mean(nPhotons(:,iSegment));

    %uncertainty value
    unc(:,iSegment) = expDelQ(:,iSegment).*expDelP(:,iSegment);
    meanUnc(iSegment) = mean(unc(:,iSegment));
    
    [~,I] = max(expQ(:,iSegment));    
    delQ(iSegment) = (expDelQ(I,iSegment));  %chose phase for evaluation of uncertainties where Q is max.
    delP(iSegment) = (expDelP(I,iSegment));    

    % Plot
    close all;
    x = meanTheta(:,iSegment);
    plot(x, expP(:,iSegment), x, expP2(:,iSegment), x, expQ(:,iSegment), x,...
        expQ2(:,iSegment), x, expDelQ(:,iSegment), x, expDelP(:,iSegment),...
        x,  unc(:,iSegment), x, nPhotons(:,iSegment), 'linewidth', 2);

    hold on;
    plot(x,Norm^2*ones(length(x)),'k-','lineWidth',0.5);
    xlabel('Phase \theta');
    axis([0 2*pi min(min(meanX))-0.5 max(max(meanX2))+0.5]);
    set(0,'DefaultLegendInterpreter','latex');
    set(0,'DefaultTextInterpreter','latex');
    legend('$<P>$', '$<P^{2}>$', '$<Q>$', '$<Q^{2}>$', '$\Delta Q$',...
        '$\Delta P$', '$\Delta Q \cdot \Delta P$ ', '$<n>$ ','uncertainty limit', 'location', 'best');
    title(strcat('Mean Photon number =','\, ',num2str(meanN(iSegment))));

    % Save plot
    print(strcat('expect-',filename,'-Segment',num2str(iSegment),'-',...
        strrep(num2str(meanN(iSegment)),'.',','),'photons'), '-dpng');

end % iSegment 

end % computeExpectations

