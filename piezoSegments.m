function [ segments ] = piezoSegments( timestamps, sourceArray )
%PIEZOSEGMENTS(timestamps,sourceArray) divides the source Arrays into
%segments corresponding to single flanks of the piezo triangular
%modulation. It looks for long gaps in the TIMESTAMPS array.

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

if nRows>1 % 2D input arrays
    segments = NaN(nRows,lSegments,nSegments);
    for iGap=1:nSegments
        seg = sourceArray(:,gapStarts(iGap)+1:gapStarts(iGap+1));
        segments(:,:,iGap) = [seg NaN(nRows, lSegments-size(seg,2))];
    end
else % 1D input arrays
    segments = NaN(lSegments,nSegments);
    for iGap=1:nSegments
        seg = sourceArray(:,gapStarts(iGap)+1:gapStarts(iGap+1));
        segments(:,iGap) = [seg NaN(nRows, lSegments-length(seg))];
    end
end

end