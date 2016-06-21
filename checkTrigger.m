function [ triggereffect ] = checkTrigger( data8bit, filename )
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

maxRandomJitter = 3;
triggereffect = zeros(maxRandomJitter,1);

for randomJitter=0:maxRandomJitter
% RANDOM JITTER
[rows,columns] = size(data8bit);
assert(rows < 67000 & columns < 180 ,'Too big data matrix for reasonable calculation time!');
rdata8bit = zeros(rows-randomJitter,columns);
for column=1:columns
    shift = randi(randomJitter+1);
    rdata8bit(:,column)=data8bit(1+(shift-1):end-(randomJitter-shift+1),1);
end
data8bit = rdata8bit;

% CALCULATE TRIGGEREFFECT
[rows,columns] = size(data8bit);

% Check width of autocorrelation peaks in single data segments
numLags = 2000;
assert(rows>numLags,'Too few data rows!');
acf=zeros(numLags,columns);
for column=1:columns
    acf(:,column) = autocorr(single(data8bit(:,column)),numLags-1);
end
sharpacf = mean(transpose(acf(:,:)));
[~,~,sharpwidths,~] = findpeaks(sharpacf,'MinPeakHeight',0);

% Check width of autocorrelation peaks after averaging all data segments
meandata = mean(transpose(data8bit(:,:)));
broadacf = autocorr(meandata,numLags-1);
[~,~,broadwidths,~] = findpeaks(broadacf,'MinPeakHeight',0);

% Compare widths
triggereffect(randomJitter+1) = mean(broadwidths)-mean(sharpwidths);
end

%CREATE PLOT
x = 0:maxRandomJitter;
bar(x, triggereffect);
xlabel('Additional uniformly distributed jitter [Samples]');
ylabel('Jitter Effect [a.u.]');

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

