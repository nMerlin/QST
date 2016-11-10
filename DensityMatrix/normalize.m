function [ normalizedMatrix ] = normalize( matrix )
%NORMALIZE normalizes a matrix to unity trace

normalizedMatrix = matrix / trace(matrix);

end

