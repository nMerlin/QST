function [ avPr ] = phaseAveragedDist( WF )
%PHASEAVERAGEDDIST Simulate phase-averaged quadrature measurement
%   
%   Arguments:
%       WF - discretized Wigner function of the target state

qRange = -20:0.125:20;

projections = zeros(360,length(qRange)); % projections of rotated WFs
for nTheta = 1:360
    % rotate and project
    rotWF = imrotate(WF,nTheta,'crop');
    projections(nTheta,:) = real(sum(rotWF));
end

avPr = sum(projections)/360; % Factor 64 has to be somewhere in WF...
plot(avPr);

end

