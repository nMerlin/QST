function [theta, ys] = computePhase(Xa,Xb,piezoSign,varargin)
%COMPUTEPHASE Reconstruct phase between Xa and Xb which can be obtained with
% prepare3ChData.
%
% Output Arguments:
%   theta - Reconstructed phase values
%   ys - smoothed crosscorrelation between Xa and Xb. This is needed for
%   segmentValidation.

%% Constants
% For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
periodsPerSeg = 1.2;

%% Handle optional input arguments
nVarargin = length(varargin);
optArgs = {'noplot'};
optArgs(1:nVarargin) = varargin;
[plotArg] = optArgs{:};

%% Reconstruct the phase for each piezo segment
ys = smoothCrossCorr(Xa,Xb);
[nPoints,nSegments] = size(ys);
ys = reshape(smooth(ys,5001,'moving'),[nPoints nSegments]);
periodLength = length(ys)*periodsPerSeg;
theta = zeros(nPoints,nSegments);
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
        title({'Mousebutton: Reject','Keyboard: Accept'});
        key = waitforbuttonpress;
        if key == 0 %mouse
            theta(:,iSeg) = NaN(nPoints,1);
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
    % closer to the boundary than 2% of a period, then it is rejected, if
    % its value is lower than 95% of the corresponding second extremum.
    % These tests work only if there are at least three peaks.
    %
    % Left boundary:
    if length(pks)>= 3
        if locs(1)<0.02*periodLength
            if (pks(1)>0 && (pks(3)-pks(1))/abs(pks(3))>0.05) || ...
                    (pks(1)<0 && (pks(1)-pks(3))/abs(pks(3))>0.05)
                locs = locs(2:end);
                pks = pks(2:end);
            end
        end
    end
    % Right boundary:
    if length(pks)>= 3
        if (length(y)-locs(end))<0.02*periodLength
            if (pks(end)>0 && (pks(end-2)-pks(end))/abs(pks(end-2))>0.05) ||...
                    (pks(end)<0 && (pks(end)-pks(end-2))/abs(pks(end-2))>0.05)
                locs = locs(1:end-1);
                pks = pks(1:end-1);
            end
        end
    end    
    nTurningPoints = length(locs);
    assert(nTurningPoints>1, 'Not enough turning points encountered!');
    
    %% Account for extrema lying directly on a boundary (1)
    % If an extremal point is the first or last point in the data set, then
    % the _findpeaks_ function won't detect it. Therefore, we have to catch
    % these exceptions here. If the boundary extremum is lower than the
    % corresponding second extremum, nothing happens, otherwise it will be
    % added to the list of peaks.
    %
    % Left boundary:
    if (pks(1)<0 && y(1)>pks(2)) || (pks(1)>0 && y(1)<pks(2))
        locs = [1;locs];
        pks = [y(1);pks];
    end
    % Right boundary:
    if (pks(end)<0 && y(end)>pks(end-1)) || ...
            (pks(end)>0 && y(end)<pks(end-1))
        locs = [locs;length(y)];
        pks = [pks;y(end)];
    end
    
    %% Account for extrema lying directly on a boundary (2)
    % Similar to the previous correction, _findpeaks_ could detect a peak,
    % but the value directly at the boundary is higher than this peak. In
    % this case, we have to replace the extremum with the boundary point.
    %
    % Left boundary:
    if (pks(1)<0 && y(1)<pks(1)) || (pks(1)>0 && y(1)>pks(1))
        pks(1) = y(1);
        locs(1) = 1;
    end
    % Right boundary:
    if (pks(end)<0 && y(end)<pks(end)) || (pks(end)>0 && y(end)>pks(end))
        pks(end) = y(end);
        locs(end) = length(y);
    end
    
    if strcmp(plotArg,'plot')
        plot(y); hold on;
        plot(locs,pks,'ro'); hold off;
        legend(num2str(iSeg));
        waitforbuttonpress;
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
            normDiff = abs(pksDiff(1));
            maxValue = max([pks(1),pks(2)]);
        elseif iPart == nTurningPoints
            range = (locs(end)):length(theta(:,iSeg));
            normDiff = abs(pksDiff(end));
            maxValue = max([pks(end),pks(end-1)]);
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
