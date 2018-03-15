function [valsConv,unitLabel] = convenientUnits(valsOrig,unitSI)
%CONVENIENTUNITS Converts SI-unit values into units with prefix (e.g. ns)
%
% Input Arguments:
%   valsOrig: Vector of input values to convert.
%   unitSI: String with original SI unit
%
% Output Arguments:
%   valsConv: Vector with converted values.
%   unitLabel: String with original SI unit and Prefix.

orderOfMagnitude = floor(log10(max(valsOrig)));
orderOfPrefix = floor(orderOfMagnitude/3);
factor = 10^(orderOfPrefix*3);
prefix = '';

switch orderOfPrefix
    case 1
        prefix = 'k';
    case 2
        prefix = 'M';
    case 3
        prefix = 'G';
    case 4
        prefix = 'T';
    case 5
        prefix = 'P';
    case -1
        prefix = 'm';
    case -2
        prefix = 'u';
    case -3
        prefix = 'n';
    case -4
        prefix = 'p';
    case -5
        prefix = 'f';
end

valsConv = valsOrig/factor;
unitLabel = [prefix,unitSI];

end

