function [ locs, v ] = pointwiseVariance( data )
%POINTWISEVARIANCE 
s = size(data);
v = zeros(s(1),1);

for i=1:s(2)
    data(:,i) = data(:,i) - mean(data(:,i));
end

for i=1:s(1)
    v(i) = var(data(i,:));
end

[pks, locs] = findpeaks(v,'MinPeakHeight',mean(v));
findpeaks(v,'MinPeakHeight',mean(v));

end

