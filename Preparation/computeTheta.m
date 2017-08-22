function [ X, theta ] = computeTheta( X, piezoSign, varargin )
%COMPUTETHETA Computes X and THETA ready for the reconstruction algorithm
%
%   X should have the size [nPulses, nRecords, nSegments]

% Optional input arguments
verbose = 0;
sinefit = 0;
quiet = 'notquiet';
if nargin > 2
    for i = 1:(nargin-2)
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

XOld = X;
[nPulses, nRecords, nSegments] = size(XOld);
[X, theta] = deal(zeros(nPulses * nRecords, nSegments));
dispstat('Computing phase values ...','timestamp',...
        'keepthis',quiet);
for iSeg = 1:nSegments
    % Data to operate on
    X(:,iSeg) = reshape(XOld(:,:,iSeg), nPulses * nRecords,1);
    
    % Average over recordings and create fit input x/y data
    yFit = mean(XOld(:,:,iSeg));
    yFit=smooth(yFit,5)';   %smoothing of yFit decreases "Dips" in theta
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
        
        yFitNew = yFit - (min(yFit) + 0.5*(max(yFit) - min(yFit))); %subtract offset for finding peaks
         
        peakHeight = 0.3 * (max(yFitNew));
        peakHeightMin = 0.3 * (max(-yFitNew));
        peakDistance = 0.3 * length(yFit);
        peakWidth = 0.01 * length(yFit);
      
%         findpeaks(yFitNew,'MinPeakHeight',peakHeight,...
%             'MinPeakDistance',peakDistance,'MinPeakWidth',peakWidth);
%         hold on;
%         findpeaks(-yFitNew,'MinPeakHeight',peakHeightMin,...
%             'MinPeakDistance',peakDistance,'MinPeakWidth',peakWidth);
        [~,maxlocs] = findpeaks(yFitNew,'MinPeakHeight',peakHeight,...
            'MinPeakDistance',peakDistance,'MinPeakWidth',peakWidth);
        [~,minlocs] = findpeaks(-yFitNew,'MinPeakHeight',peakHeightMin,...
            'MinPeakDistance',peakDistance,'MinPeakWidth',peakWidth);
        
        maxpks = yFit(maxlocs);
        minpks = yFit(minlocs);
        
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
        [minVal,I] = min(yFit(1:locs(1)));
        if pks(1)<0 && minVal<pks(1)
            locs(1) = I;
            pks(1) = minVal;
        end
        % First peak is a maximum and left extremum is a higher maximum
        [maxVal,I] = max(yFit(1:locs(1)));
        if pks(1)>0 && maxVal>pks(1)
            locs(1) = I;
            pks(1) = maxVal;
        end
        % First peak is a minimum and left extremum is a maximum
        if pks(1)<0
            leftex = max(yFit(1:locs(1)));
        % First peak is a maximum and left extremum is a minimum
        else
            leftex = min(yFit(1:locs(1)));
        end
        
        %% Look for extrema on right boundary
        % Last peak is a minimum and right extremum is a lower minimum
        [minVal,I] = min(yFit(locs(end):end));
        if pks(end)<0 && minVal<pks(end)
            locs(end) = I;
            pks(end) = minVal;
        end
        % Last peak is a maximum and right extremum is a higher maximum
        [maxVal,I] = max(yFit(locs(end):end));
        if pks(end)>0 && maxVal>pks(end)
            locs(end) = I;
            pks(end) = maxVal;
        end
        % Last peak is a minimum and right extremum is a maximum
        if pks(end)<0
            rightex = max(yFit(locs(end):end));
        % Last peak is a maximum and right extremum is a minimum
        else
            rightex = min(yFit(locs(end):end));
        end

        %% Loop over all visible flanks
        smallTheta = zeros(length(yFit),1);
        % direction of the first visible flank; also account for the
        % different directions of the piezo movement
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
            
            % Correct X-Offset
            offset = (2*maxValue-normDiff)/2;
            Xrange = ((range(1)-1)*nPulses+1):(range(end)-1*nPulses);
            X(Xrange) = X(Xrange) - offset;
            
            % Scale y-Values to interval [-1;1] for asin
            y = yFit(range);
            ynorm = 2*(y-maxValue)/normDiff + 1;
            
            % Correct for machine precision
            [~,iMax] = max(ynorm);
            ynorm(iMax) = ynorm(iMax) - 2*eps;
            [~,iMin] = min(ynorm);
            ynorm(iMin) = ynorm(iMin) + 2*eps;
            
            
%             hold on;
%             plot(ynorm);
%             hold on;
%             
            % Calculate phases
            if s==1
                smallTheta(range) = asin(ynorm);
            else
                smallTheta(range) = pi - asin(ynorm);
            end
                    
            if ( piezoSign == 1)
                if ( ss == 1)
                     smallTheta(range) = smallTheta(range)+2*pi*floor(iPart/2);                
                else
                     smallTheta(range) = smallTheta(range)+...
                     2*pi*floor((iPart+1)/2);
                end
            else
                if ( ss == 1 )
                     smallTheta(range) = smallTheta(range)-2*pi*floor((iPart+1)/2);                
                else
                     smallTheta(range) = smallTheta(range)-...
                     2*pi*floor(iPart/2);
                end
                
            end
            s = s * (-1);
        end
        
        
%         plot(smallTheta);
%         hold on;
        
        assert(isreal(smallTheta),...
            ['Not all phase values are real in Segment ' num2str(iSeg) '.']);
        
        % Calculate phase values from inperpolated "smallTheta"
        xSample = 1 : nPulses * nRecords;
        theta(:,iSeg) = mod(interp1(xFit(~isnan(xFit)),...
            smallTheta(~isnan(smallTheta)),xSample,'spline','extrap'),2*pi);
        theta(isnan(X(:,iSeg)),iSeg) = NaN;

        %subtract offset
        [~,I] = sort(theta(:,iSeg));
        Xsorted = X(I,iSeg);
        XsortedSmoothed = smooth(Xsorted,20,'moving');
        [~,Imax] = max(XsortedSmoothed);
        [~,Imin] = min(XsortedSmoothed);
        span = 100;
        if Imin<Imax
            minValue = mean(Xsorted(Imin:Imin+span));
            maxValue = mean(Xsorted(Imax-span:Imax));
        else
            minValue = mean(Xsorted(Imin-span:Imin));
            maxValue = mean(Xsorted(Imax:Imax+span));
        end
        offset = minValue+0.5*(maxValue-minValue);
        X(:,iSeg)=X(:,iSeg)-offset;

    end
end

%Concatenate the segments 
% theta = reshape(theta,size(theta,1)*size(theta,2) ,1);
% X = reshape(X,size(X,1)*size(X,2) ,1);

end

