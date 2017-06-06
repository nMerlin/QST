function [ avPr ] = phaseAveragedDist( WF )
%PHASEAVERAGEDDIST Simulate phase-averaged quadrature measurement
%   
%   Arguments:
%       WF - discretized Wigner function of the target state

qRange = -20:0.125:20;

% Convert WF into (p,q,z) dataset called BIGWF
bigWF = zeros(3,length(WF)^2);
for p = 1:length(qRange)
    for q = 1:length(qRange)
        bigWF(1,(p-1)*length(qRange)+q) = qRange(p);
        bigWF(2,(p-1)*length(qRange)+q) = qRange(q);
        bigWF(3,(p-1)*length(qRange)+q) = WF(p,q);
    end
end

projections = zeros(360,length(qRange)); % projections of rotated WFs
for nTheta = 1:360
    % rotate and project
    rotWF = imrotate(WF,nTheta,'crop');
    projections(nTheta,:) = real(sum(rotWF));
end

avPr = sum(projections)/(360*64); % Factor 64 has to be somewhere in WF...
plot(avPr);

end

