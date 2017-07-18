function [theta, selSeg] = computePhase(ys,piezoSign,varargin)
%COMPUTEPHASE Reconstruct phase from smoothed cross-correlation data
%
% Input Arguments:
%   ys = smoothCrossCorr(Xa,Xb,varargin)
%
% Output Arguments:
%   theta - Reconstructed phase values
%   selSeg - the reconstruction was successfull for these segments

%% Global parameters
% For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
periodsPerSeg = 1.2;
periodLength = length(ys)*periodsPerSeg;

%% Handle optional input arguments
nVarargin = length(varargin);
optArgs = {'noplot'};
optArgs(1:nVarargin) = varargin;
[plotArg] = optArgs{:};

[nPoints,nSegments] = size(ys);
theta = zeros(nPoints,nSegments);
selSeg = ones(1,nSegments);
for iSeg = 1:nSegments
    y = ys(:,iSeg);
    %% Parameters for the function _findpeaks_
    % To reconstruct the phase, the algorithm has to identify maxima and
    % minima in the data. This is done by the function _findpeaks_.
    % However, the parameters have to be chosen carefully.
    peakOpts.MinPeakDistance = 0.6 * length(y)/periodsPerSeg;
    % _MinPeakDistance_ is by far the most important parameter. It
    % determines how far away of each other the found peaks must be. In our
    % case, the data should exhibit a certain periodicity and we want to
    % know where the maximum and minimum in each period is. Therefore, a
    % _MinPeakDistance_ of roughly 50% of the period should do the trick.
    peakOpts.MinPeakHeight = 0.5 * max(y);
    % Because of different noise sources and instabilities, the distance
    % between maxima and minima in _y_ is uncertain to some degree.
    % _MinPeakHeight_ ensures, that only peaks above a certain threshold
    % will be found. Because between two maxima there should be a minimum,
    % _MinPeakHeight_ mitigates the effects of too long periods.
    
    %% Find peaks for flank recognition
    if strcmp(plotArg,'plot')
        findpeaks(y,peakOpts,'Annotate','extents'); hold on;
        findpeaks(-y,peakOpts,'Annotate','extents'); hold off;
        key = waitforbuttonpress;
        if key == 0
            theta(:,iSeg) = NaN(nPoints,1);
            selSeg(iSeg) = false;
            continue
        end
    end
    [~,maxlocs] = findpeaks(y,peakOpts);
    [~,minlocs] = findpeaks(-y,peakOpts);
    maxpks = y(maxlocs);
    minpks = y(minlocs);
    assert(abs(length(maxpks)-length(minpks))<2,...
            strcat('Too many maxima or minima detected in Segment', ...
            num2str(iSeg),'!'));
    
    %% Sort peaks (assumption: we only see "global" maxima and minima)
    [locs, I] = sort([maxlocs;minlocs]);
    pks = [maxpks;minpks];
    pks = pks(I);
    
    %% Account for wrongly detected peaks close to the boundaries
    % An extremum very close to the boundary could be the result of local
    % fluctuations instead of the piezo modulation. If such an extremum is
    % closer to the boundary than 5% of a period, then it is rejected, if
    % its value is not within a 5% range of the second extremum in the
    % direction of 0.
    if locs(1)<0.05*periodLength && abs((pks(1)-pks(3))/pks(3))>0.05
        if (pks(1)>0 && (pks(2)-pks(1))>0.5) || ...
                (pks(1)<0 && (pks(1)-pks(2))>0.5)
            locs = locs(2:end);
            pks = pks(2:end);
        end
    end
    
    nTurningPoints = length(locs);
    assert(nTurningPoints>1, 'Not enough turning points encountered!');
    
    %% Catch ...
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
    
    plot(y); hold on;
    plot(locs,pks,'ro'); hold off;
    waitforbuttonpress;
    
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
            range = (locs(end)):length(theta(:,iSeg));
            normDiff = max(abs(pksDiff(end)),abs(rightex-pks(end)));
            maxValue = max([pks(end),pks(end-1),rightex]);
        else
            range = (locs(iPart)):(locs(iPart+1));
            normDiff = abs(pksDiff(iPart));
            maxValue = max(pks(iPart),pks(iPart+1));
        end

        % Scale y-Values to interval [-1;1] for asin
        ynorm = y(range);
        ynorm = 2*(ynorm-maxValue)/normDiff + 1;

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
    
    if ~isreal(theta(1:end,iSeg))
        theta(:,iSeg) = NaN(nPoints,1);
        selSeg(iSeg) = false;
        continue
    end
    assert(isreal(theta),...
        ['Not all phase values are real in Segment ' num2str(iSeg) '.']);
end % iSeg

theta = mod(theta,2*pi);

end % function
