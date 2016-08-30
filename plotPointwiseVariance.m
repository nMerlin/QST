function [  ] = plotPointwiseVariance( data, filename )
% PLOTPOINTWISEVARIANCE calculates the pointwise variance of DATA, plots
% the pointwise variance data and highlights the maxima that are above the
% mean value. The columns of DATA are the datasets.
%
%   PLOTPOINTWISEVARIANCE(DATA) Plots the pointwise variance of DATA.
%
%   PLOTPOINTWISEVARIANCE(DATA, FILENAME) Plots the pointwise variance of
%   DATA and outputs the resulting graph to FILENAME. Does the filename
%   include an extenstion, this will be used as the filetype. If not, 'png'
%   is the standard filetype.
%   
%   See also: POINTWISEVARIANCE

[~,pvar] = pointwiseVariance(data);

assert(isvector(pvar),'PVAR is not a vector!');
findpeaks(pvar,'MinPeakHeight',mean(pvar),'MinPeakDistance',10);

switch nargin
    case 2
        assert(ischar(filename),'FILENAME is not a character array!');
        [~,~,extension]=fileparts(filename);
        if not(isempty(extension))
            extension = strsplit(extension,'.');
            filetype = char(extension(2));
        else
            filetype = 'png';
        end
        
        if strcmp(filetype,'jpg')
            filetype = 'jpeg';
        end
        filetype = ['-d' filetype];
        
        assert(ischar(filetype),'FILETYPE is not a character array!');
        print(filename,filetype);
end

end