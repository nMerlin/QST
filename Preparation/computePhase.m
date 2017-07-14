function theta = computePhase(ys,piezoSign,varargin)
%COMPUTEPHASE Reconstruct phase from smoothed cross-correlation ys
%
% Input Arguments:
%   ys = smoothCrossCorr(Xa,Xb,varargin)

%% Handle optional input arguments
nVarargin = length(varargin);
optArgs = {'noplot'};
optArgs(1:nVarargin) = varargin;
[plotArg] = optArgs{:};

%% Parameters for the function _findpeaks_
peakOpts.MinPeakHeight = 0.3 * max(ys);
peakOpts.MinPeakDistance = 0.3 * length(ys);
peakOpts.MinPeakWidth = 0.01 * length(ys);

[nPoints,nSegments] = size(ys);
theta = zeros(nPoints,nSegments);

for iSeg = 1:nSegments
    y = ys(:,iSeg);
    
    %% Find peaks for flank recognition
    if strcmp(plotArg,'plot')
        findpeaks(y,peakOpts); hold on;
        findpeaks(-y,peakOpts); hold off;
    end
    [~,maxlocs] = findpeaks(y,peakOpts);
    [~,minlocs] = findpeaks(-y,peakOpts);
    maxpks = y(maxlocs);
    minpks = y(minlocs);
    assert(abs(length(maxpks)-length(minpks))<2,...
            strcat('Too many maxima or minima detected in Segment', ...
            num2str(iSeg),'!'));
    
    %% Sort peaks (assumption: we only see "global" maxima and minima)
    [locs, I] = sort([maxlocs minlocs]);
    nTurningPoints = length(locs);
    assert(nTurningPoints>1, 'Not enough turning points encountered!');
    pks = [maxpks minpks];
    pks = pks(I);
    
    %% Look for extrema on left boundary
    % First peak is a minimum and left extremum is a lower minimum
    [minVal,I] = min(y(1:locs(1)));
    if pks(1)<0 && minVal<pks(1)
        locs(1) = I;
        pks(1) = minVal;
    end
    % First peak is a maximum and left extremum is a higher maximum
    [maxVal,I] = max(y(1:locs(1)));
    if pks(1)>0 && maxVal>pks(1)
        locs(1) = I;
        pks(1) = maxVal;
    end
    % First peak is a minimum and left extremum is a maximum
    if pks(1)<0
        leftex = max(y(1:locs(1)));
    % First peak is a maximum and left extremum is a minimum
    else
        leftex = min(y(1:locs(1)));
    end
        
    %% Look for extrema on right boundary
    % Last peak is a minimum and right extremum is a lower minimum
    [minVal,I] = min(y(locs(end):end));
    if pks(end)<0 && minVal<pks(end)
        locs(end) = I;
        pks(end) = minVal;
    end
    % Last peak is a maximum and right extremum is a higher maximum
    [maxVal,I] = max(y(locs(end):end));
    if pks(end)>0 && maxVal>pks(end)
        locs(end) = I;
        pks(end) = maxVal;
    end
    % Last peak is a minimum and right extremum is a maximum
    if pks(end)<0
        rightex = max(y(locs(end):end));
    % Last peak is a maximum and right extremum is a minimum
    else
        rightex = min(y(locs(end):end));
    end
    
    %% Loop over all visible flanks
    % _ss_ accounts for the direction of the first visible flank and for
    % the different directions of the piezo movement
    pksDiff = -diff(pks);
    ss = sign(pksDiff(1))*piezoSign;
    s = ss;
    for iPart = 0:nTurningPoints
        % Normalize to interval [-1;1]
        if iPart == 0
            range = 1:locs(1);
            normDiff = max(abs(pksDiff(1)),abs(leftex-pks(1)));
            maxValue = max([pks(1),pks(2),leftex]);
        elseif iPart == nTurningPoints
            range = (locs(end)):length(smallTheta);
            normDiff = max(abs(pksDiff(end)),abs(rightex-pks(end)));
            maxValue = max([pks(end),pks(end-1),rightex]);
        else
            range = (locs(iPart)):(locs(iPart+1));
            normDiff = abs(pksDiff(iPart));
            maxValue = max(pks(iPart),pks(iPart+1));
        end

        % Scale y-Values to interval [-1;1] for asin
        y = yFit(range);
        ynorm = 2*(y-maxValue)/normDiff + 1;

        % Correct for machine precision
        [~,iMax] = max(ynorm);
        ynorm(iMax) = ynorm(iMax) - 2*eps;
        [~,iMin] = min(ynorm);
        ynorm(iMin) = ynorm(iMin) + 2*eps;
           
        % Calculate phases
        if s==1
            theta(range,iSeg) = asin(ynorm);
        else
            theta(range,iSeg) = pi - asin(ynorm);
        end

        if ( piezoSign == 1)
            if ( ss == 1)
                 theta(range,iSeg) = theta(range,iSeg)+2*pi*floor(iPart/2);                
            else
                 theta(range,iSeg) = theta(range,iSeg)+ ...
                     2*pi*floor((iPart+1)/2);
            end
        else
            if ( ss == 1 )
                 theta(range,iSeg) = theta(range,iSeg)- ...
                     2*pi*floor((iPart+1)/2);                
            else
                 theta(range,iSeg) = theta(range,iSeg)- ...
                 2*pi*floor(iPart/2);
            end
        end
        s = s * (-1);
    end % iPart
    
    assert(isreal(theta),...
        ['Not all phase values are real in Segment ' num2str(iSeg) '.']);
end % iSeg

end % function
