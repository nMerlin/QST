function [ randomNumbers ] = distributeRandomly( n, x, probabilityDensity )
%DISTRIBUTERANDOMLY Compute random numbers for an arbitrary distribution
%
%   N: amount of random numbers to generate
%   X: x-axis of the probability density function
%   PROBABILITYDENSITY: values of the probability density function
%
%   Adapted from:
%   http://matlabtricks.com/
%   post-44/generate-random-numbers-with-a-given-distribution

% normalize the probability density function
probabilityDensity = probabilityDensity / sum(probabilityDensity);

% compute cumulative distribution function
cumulativeDistribution = cumsum(probabilityDensity);

% remove non-unique elements
[cumulativeDistribution, mask] = unique(cumulativeDistribution);
x = x(mask);

% create an array of N random numbers
randomNumbers = rand(1,n);

% inverse interpolation to compute P(x) -> x projection
randomNumbers = interp1(cumulativeDistribution, x, randomNumbers);

end

