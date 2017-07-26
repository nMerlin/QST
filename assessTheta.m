function [variance, relVariance] = assessTheta(theta)


%% Create and plot histogram
% The number of bins is 1000 by default.
maxVal = max(-min(theta),max(theta));
xHist = linspace(-maxVal,maxVal,1000);
disc = min(diff(xHist));
histEdges = (-maxVal-disc/2):disc:(maxVal+disc/2);


h = histogram(theta,histEdges,'Normalization','pdf');
%plotting
h.EdgeColor = 'b';

axis([0 2*pi 0 max(h.Values)]);
xlabel('\theta');
ylabel('probability density');

%values
variance = var(h.Values);

relVariance = variance/mean(h.Values);
end
