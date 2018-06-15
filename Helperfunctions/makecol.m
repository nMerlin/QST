function vector = makecol(vector)
%MAKECOL converts any input vector into a column vector

if isrow(vector)
    vector = vector';
end

end

