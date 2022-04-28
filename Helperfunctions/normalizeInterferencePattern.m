function [ynew,locs,maxlocs,minlocs] = normalizeInterferencePattern(y,numberOfPeriods)
%% normalize the interference pattern to absolute maximum and minimum 
periodLength = length(y)/numberOfPeriods;
peakOptsMax.MinPeakDistance = periodLength;
peakOptsMin.MinPeakDistance = periodLength;
[~,maxlocs] = findpeaks(y, peakOptsMax);
[~,minlocs] = findpeaks(-y,peakOptsMin);
maxpks = y(maxlocs);
minpks = y(minlocs);
%sort peaks
[locs, I] = sort([maxlocs;minlocs]);
pks = [maxpks;minpks];
pks = pks(I);
absMax = mean(maxpks);
absMin = mean(minpks);

% %% Account for wrongly detected peaks close to the boundaries
%     % An extremum very close to the boundary could be the result of local
%     % fluctuations instead of the piezo modulation. If such an extremum is
%     % closer to the boundary than 2% of a period, then it is rejected, if
%     % its value is lower than 95% of the corresponding second extremum.
%     % These tests work only if there are at least three peaks.
%     %
%     % Left boundary:
%     if length(pks)>= 3
%         if locs(1)<0.02*periodLength
%             if (pks(1)>0 && (pks(3)-pks(1))/abs(pks(3))>0.05) || ...
%                     (pks(1)<0 && (pks(1)-pks(3))/abs(pks(3))>0.05)
%                 locs = locs(2:end);
%                 pks = pks(2:end);
%             end
%         end
%     end
%     % Right boundary:
%     if length(pks)>= 3
%         if (length(y)-locs(end))<0.02*periodLength
%             if (pks(end)>0 && ...
%                 (pks(end-2)-pks(end))/abs(pks(end-2))>0.05) || ...
%                 (pks(end)<0 && (pks(end)-pks(end-2))/abs(pks(end-2))>0.05)
%                 locs = locs(1:end-1);
%                 pks = pks(1:end-1);
%             end
%         end
%     end    
%     
% %     assert(nTurningPoints>1, ...
% %         ['Not enough turning points encountered in Segment ', ...
% %         num2str(iSeg),'!']);
%     
    %% Account for extrema lying directly on a boundary (1)
    % If an extremal point is the first or last point in the data set, then
    % the _findpeaks_ function won't detect it. Therefore, we have to catch
    % these exceptions here. If the boundary extremum is lower than the
    % corresponding second extremum, nothing happens, otherwise it will be
    % added to the list of peaks.
    %
    % Left boundary:
%     if (pks(1)<0 && y(1)>pks(2)) || (pks(1)>0 && y(1)<pks(2))
%         locs = [1;locs];
%         pks = [y(1);pks];
%     end
%     % Right boundary:
%     if (pks(end)<0 && y(end)>pks(end-1)) || ...
%             (pks(end)>0 && y(end)<pks(end-1))
%         locs = [locs;length(y)];
%         pks = [pks;y(end)];
%     end
    
    %% Account for extrema lying directly on a boundary (2)
    % Similar to the previous correction, _findpeaks_ could detect a peak,
    % but the value directly at the boundary is higher than this peak. In
    % this case, we have to replace the extremum with the boundary point.
    %
%     % Left boundary:
%     if (pks(1)<0 && y(1)<pks(1)) || (pks(1)>0 && y(1)>pks(1))
%         pks(1) = y(1);
%         locs(1) = 1;
%     end
%     % Right boundary:
%     if (pks(end)<0 && y(end)<pks(end)) || (pks(end)>0 && y(end)>pks(end))
%         pks(end) = y(end);
%         locs(end) = length(y);
%     end
    

     %% Loop over all visible flanks

pksDiff = -diff(pks);
nTurningPoints = length(locs);
ynew = y;
for iPart = 0:nTurningPoints
    % Normalize to interval [-1;1]
    if iPart == 0
        range = 1:locs(1);
        normDiff = abs(pksDiff(1));
        maxValue = max([pks(1),pks(2)]);
    elseif iPart == nTurningPoints
        range = (locs(end)):length(y);
        normDiff = abs(pksDiff(end));
        maxValue = max([pks(end),pks(end-1)]);
    else
        range = (locs(iPart)):(locs(iPart+1));
        normDiff = abs(pksDiff(iPart));
        maxValue = max(pks(iPart),pks(iPart+1));
    end

    % Scale y-Values to interval [absMax;absMin] 
    ynorm = y(range);
    ynorm = (absMax-absMin)*(ynorm-maxValue)/normDiff + absMax;

    % Correct for machine precision
    [~,iMax] = max(ynorm);
    ynorm(iMax) = ynorm(iMax) - 2*eps;
    [~,iMin] = min(ynorm);
    ynorm(iMin) = ynorm(iMin) + 2*eps;
    ynew(range) = ynorm;
end
end