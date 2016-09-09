function [  ] = plotCorrelationValues( data8bit, windowSize, filename )
%PLOTCORRELATIONVALUES Summary of this function goes here
%   Detailed explanation goes here

% Parameters & Variables
jRange = 0:10;
correlations = zeros(1,length(jRange));

% Capture integration centers
[locs,~] = pointwiseVariance(data8bit);

parfor j=1:length(jRange)
    correlations(j) = correlation(jRange(j),data8bit,locs,windowSize);
end

%%% Plotting
close all;
bar(jRange, correlations);
for i1=1:numel(correlations)
    text(jRange(i1),correlations(i1),num2str(correlations(i1),'%0.2f'), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top');
end
ylabel('Correlation of Q_i and Q_{i+j}');
xlabel('j');

% Saving
switch nargin
    case 3
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

