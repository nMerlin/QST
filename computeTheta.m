function [ X, theta ] = computeTheta( X, varargin )
%COMPUTETHETA Computes X and THETA ready for the reconstruction algorithm
%
%   X should have the size [nPulses, nRecords, nSegments]

% Optional input arguments
verbose = 0;
sinefit = 0;
quiet = 'notquiet';
if nargin > 1
    for i = 1:(nargin-1)
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

XOld = X;
[nPulses, nRecords, nSegments] = size(XOld);
[X, theta] = deal(zeros(nPulses * nRecords, nSegments));
dispstat('Calculating phase values...','timestamp',...
        'keepthis',quiet);
for iSeg = 1:nSegments
    % Data to operate on
    X(:,iSeg) = reshape(XOld(:,:,iSeg), nPulses * nRecords,1);
    
    % Average over recordings and create fit input x/y data
    yFit = mean(XOld(:,:,iSeg));
    xFit = (1:nRecords) * nPulses - round(nPulses/2);
    xFit(isnan(yFit)) = NaN;
    
    % You can choose between different phase reconstruction schemes
    if sinefit == 1
        % Somehow the fit only works for the correct x magnitude
        xFitMagnitude = ceil(log10(max(xFit)));
        xFit = xFit / 10^(xFitMagnitude); % scale x-axis for fitting routine

        % Fitting
        [fitParams, ~] = fitSinusoidal(xFit, yFit, 'rmLin', 'show');

        % Correcting X for offset and linear trend
        X(:,iSeg) = X(:,iSeg) - fitParams(4) - ...
            fitParams(5) * (1 : length(X(:,iSeg)))' / 10^(xFitMagnitude);

        % Calculate phase values from fit results
        xTheta = 1 : nPulses * nRecords;
        theta(:,iSeg) = ...
            mod(2 * pi / fitParams(2) * xTheta / 10^(xFitMagnitude) + ...
            2 * pi / fitParams(3) + pi / 2, 2 * pi);
        theta(isnan(X)) = NaN;
    else
        % Method: Inverse sine function on normalized y-values
        peakHeightFactor = 0.7;
        peakDistance = 0.1 * length(yFit);
        peakWidth = 0.01 * length(yFit);
        findpeaks(yFit,'MinPeakHeight',...
            peakHeightFactor*max(yFit),'MinPeakDistance',peakDistance,...
            'MinPeakWidth',peakWidth);
        [maxpks,maxlocs] = findpeaks(yFit,'MinPeakHeight',...
            peakHeightFactor*max(yFit),'MinPeakDistance',peakDistance,...
            'MinPeakWidth',peakWidth);
        findpeaks(-yFit,'MinPeakHeight',...
            peakHeightFactor*max(-yFit),'MinPeakDistance',peakDistance,...
            'MinPeakWidth',peakWidth);
        [minpks,minlocs] = findpeaks(-yFit,'MinPeakHeight',...
            peakHeightFactor*max(-yFit),'MinPeakDistance',peakDistance,...
            'MinPeakWidth',peakWidth);
        assert(abs(length(maxpks)-length(minpks))<2,...
            'Too many maxima or minima detected!');
        
        % Sort peaks (assumption: we only see "global" maxima and minima)
        [locs, I] = sort([maxlocs minlocs]);
        nTurningPoints = length(locs);
        assert(nTurningPoints>1, 'Not enough turning points encountered!');
        
        pks = [maxpks -minpks];
        pks = pks(I);
        pksDiff = -diff(pks);
        
        % Loop over all visible flanks
        smallTheta = zeros(length(yFit),1);
        ss = sign(pksDiff(1)); % direction of the first visible flank
        s = ss;
        for iPart = 0:nTurningPoints
            % Normalize to interval [-1;1]
            if iPart == 0 % left border to first peak
                range = 1:locs(1);
                if (max(yFit(1:10))<max(pks(1),pks(2)) &&...
                        min(yFit(1:10))>min(pks(1),pks(2)))
                    normDiff = abs(pksDiff(1));
                    maxValue = max(pks(1),pks(2));
                else
                    maxValue = max(max(yFit(1:10)),pks(1));
                    normDiff = abs(maxValue-min(min(yFit(1:10)),pks(1)));
                end
            elseif iPart == nTurningPoints % last peak to right border
                range = (locs(end)+1):length(smallTheta);
                if (max(yFit((end-10):end))<max(pks(end),pks(end-1)) && ...
                        min(yFit((end-10):end))>min(pks(end),pks(end-1)))
                    normDiff = abs(pksDiff(end));
                    maxValue = max(pks(end),pks(end-1));
                else
                    maxValue = max(max(yFit((end-10):end)),pks(end));
                    normDiff = abs(maxValue-min(min(yFit((end-10):end)),...
                        pks(end)));
                end
            else
                range = (locs(iPart)+1):(locs(iPart+1));
                normDiff = abs(pksDiff(iPart));
                maxValue = max(pks(iPart),pks(iPart+1));
            end
            
            % Correct X-Offset
            offset = (2*maxValue-normDiff)/2;
            Xrange = ((range(1)-1)*nPulses+1):(range(end)-1*nPulses);
            X(Xrange) = X(Xrange) - offset;
            
            % Scale y-Values to interval [-1;1] for asin
            y = yFit(range);
            ynorm = 2*(y-maxValue)/normDiff + 1;
            
            % Correct for machine precision
            if (ynorm(end)==max(ynorm))
                ynorm(end) = ynorm(end) - 2*eps;
                ynorm(1) = ynorm(1) + 2*eps;
            else
                ynorm(end) = ynorm(end) + 2*eps;
                ynorm(1) = ynorm(1) - 2*eps;
            end
            
            % Calculate phases
            if s==1
                smallTheta(range) = asin(ynorm);
            else
                smallTheta(range) = pi - asin(ynorm);
            end
            if (ss == 1)
                smallTheta(range) = smallTheta(range)+2*pi*floor(iPart/2);
            else
                smallTheta(range) = smallTheta(range)+...
                    2*pi*floor((iPart+1)/2);
            end
            s = s * (-1);
        end
        
        assert(isreal(smallTheta),'Not all phase values are real.');
        
        % Calculate phase values from inperpolated "smallTheta"
        xSample = 1 : nPulses * nRecords;
        theta(:,iSeg) = mod(interp1(xFit(~isnan(xFit)),...
            smallTheta(~isnan(smallTheta)),xSample,'spline','extrap'),2*pi);
        theta(isnan(X(:,iSeg)),iSeg) = NaN;
    end
end

end

