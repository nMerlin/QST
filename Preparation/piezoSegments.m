function [segments, piezoSign] = piezoSegments(timestamps, sourceArray, varargin)
%PIEZOSEGMENTS(timestamps,sourceArray) divides the source Arrays into
%segments corresponding to single flanks of the piezo triangular
%modulation. It looks for long gaps in the TIMESTAMPS array.

cut = 0;
if nargin > 2
    for i = 3:nargin
        eval([varargin{i-2} '=1;']);
    end
end

% Check data format
[nRows, nColumns] = size(sourceArray);
if nColumns~=length(timestamps)
    sourceArray = sourceArray';
    nRows = nColumns;
end
assert(size(sourceArray,2)==length(timestamps),'SourceArray Mismatch!');

% Segment data
t = diff(timestamps);
gapStarts = find((max(t)-min(t))/2<t);
nSegments = length(gapStarts)-1;
lSegments = int32(max(diff(gapStarts)));

% For phase reconstruction, the direction of the first segment is important
gaps = t(gapStarts);
if gaps(1)<0.5*(max(gaps)+min(gaps))
    piezoSign = 1; % Piezo moves in the first segment from 0 to 2 um
else
    piezoSign = -1; % Piezo moves in the first segment from 2 to 0 um
end
    
if nRows>1 % 2D input arrays
    segments = NaN(nRows,lSegments,nSegments);
    for iGap=1:nSegments
        seg = sourceArray(:,gapStarts(iGap)+1:gapStarts(iGap+1));
        segments(:,:,iGap) = [seg NaN(nRows, lSegments-size(seg,2))];
    end
    % Optional: Cut piezo segments to equal length without NaN
    if cut == 1
        nDelete = max(sum(sum(isnan(segments))))/nRows;
        segments = segments(:,1:end-nDelete,:);
    end
else % 1D input arrays
    segments = NaN(lSegments,nSegments);
    for iGap=1:nSegments
        seg = sourceArray(:,gapStarts(iGap)+1:gapStarts(iGap+1));
        segments(:,iGap) = [seg NaN(nRows, lSegments-length(seg))];
    end
end

%comment: these would be the "real times" of the starting points:
% data8bitCh = data8bit(:,:,iCh);
% timevector = (1:length(data8bitCh(:)))/configSIG.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL*1e6;
% [seglength,~] = size(data8bitCh);
% starttimes = timevector(gapStarts*seglength);


end