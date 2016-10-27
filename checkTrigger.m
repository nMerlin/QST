function [ triggerEffect ] = checkTrigger( data8bit, filename )
%CHECKTRIGGER Check jitter in trigger timing of DATA8BIT.
%
%   TRIGGEREFFECT=CHECKTRIGGER(DATA8BIT) First, the autocorrelation for
%   every data segment (columns) is calculated and averaged. Then the mean
%   width of the laser peaks in this data is calculated. This number is
%   compared to the mean laser peak width in the autocorrelation function
%   of the average segment (average over all segements). TRIGGEREFFECT
%   consists of the differences between the former "sharp" width and the
%   latter "broad" width for different additional random jitters.
%
%   TRIGGEREFFECT=CHECKTRIGGER(DATA8BIT,FILENAME) Additionally to the above
%   functionality, the resulting bar chart of TRIGGEREFFECT is saved to
%   FILENAME.

NMAX_RANDOM_JITTER_SAMPLES = 3;

triggerEffect = zeros(NMAX_RANDOM_JITTER_SAMPLES,1);

% Calculate jitter values for artificial jitter
for iRandomJitter=0:NMAX_RANDOM_JITTER_SAMPLES
[nRows,nColumns] = size(data8bit);

assert(nRows < 67000 & nColumns < 180 ,'Too big data matrix for reasonable calculation time!');
jitteredData8bit = zeros(nRows-iRandomJitter,nColumns);
for iColumns=1:nColumns
    shift = randi(iRandomJitter+1);
    if iRandomJitter > 0
        jitteredData8bit(:,iColumns)=data8bit(1+(shift-1):end-(iRandomJitter-shift+1),1);
    else
        jitteredData8bit(:,iColumns)=data8bit(1+(shift-1):end-(iRandomJitter-shift+1),iColumns);
    end
end
data8bit = jitteredData8bit;

% CALCULATE TRIGGEREFFECT
[nRows,nColumns] = size(data8bit);

% Check width of autocorrelation peaks in single data segments
nLags = 2000;
assert(nRows>nLags,'Too few data rows!');
autoCorrelationMatrix=zeros(nLags,nColumns);
for iColumns=1:nColumns
    if iRandomJitter==0
        autoCorrelationMatrix(:,iColumns) = autocorr(single(data8bit(:,iColumns)),nLags-1);
    else
        autoCorrelationMatrix(:,iColumns) = autocorr(single(data8bit(:,1)),nLags-1);
    end
end
% Check width of autocorrelation peaks for every single data segment
sharpAutoCorrelationMatrix = mean(transpose(autoCorrelationMatrix(:,:)));
[~,~,sharpwidths,~] = findpeaks(sharpAutoCorrelationMatrix,'MinPeakHeight',0);

% Check width of autocorrelation peaks after averaging all data segments
meandata = mean(transpose(data8bit(:,:)));
broadacf = autocorr(meandata,nLags-1);
[~,~,broadwidths,~] = findpeaks(broadacf,'MinPeakHeight',0);

% Compare widths
triggerEffect(iRandomJitter+1) = mean(broadwidths)-mean(sharpwidths);
end

%CREATE PLOT
xValues = 0:NMAX_RANDOM_JITTER_SAMPLES;
bar(xValues(1), triggerEffect(1),'r');
hold on;
bar(xValues(2:end),triggerEffect(2:end));
xlabel('First bar: measured value; Further bars: Simulated for uniformly distributed jitter');
ylabel('Jitter Effect [a.u.]');
hold off;

%Save plot
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

