function [  ] = plotStackedWaveforms( data, filename )
% PLOTSTACKEDWAVEFORMS Plots the waveforms given in the columns of DATA on top of
% each other. Before plotting, the mean value of each waveform is
% substracted from it.
%   
%   PLOTSTACKEDWAVEFORMS(data) Plots the waveforms on top of each other.
%   
%   PLOTSTACKEDWAVEFORMS(data,filename) Plots the waveforms in data and outputs
%   the resulting graph to FILENAME. Does the filename include an
%   extenstion, this will be used as the filetype. If not, 'png' is the
%   standard filetype.

assert(ismatrix(data),'DATA is not a matrix!');
[rows, columns] = size(data);

for column=1:columns
    y = data(:,column);
    y = y - mean(y);
    plot(y,'b');
    hold on;
end

xlabel('time');
ylabel('voltage');

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