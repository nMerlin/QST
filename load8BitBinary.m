function [ data8bit, config, timestamps ] = load8BitBinary( filename )
%LOAD8BITBINARY Loads 8bit binary datafiles, the configuration file and the
%timestamp file for a single multiple recording measurement with a Spectrum
%data acquisition card.
%
%   DATA8BIT = LOAD8BITBINARY(filename) Depending on the number of
%   channels, DATA8BIT is a 2D- or 3D- array with dimensions [columns,
%   rows, channels] that contains a single measured segment of one channel
%   in a column. To extract the number of channels, a configuration file
%   filename.cfg is necessary.
%
%   [DATA8BIT, CONFIG, TIMESTAMPS] = LOAD8BITBINARY(filename) Additionally
%   to the previously discussed array DATA8BIT, the structure CONFIG
%   consists of the configuration data and TIMESTAMPS is a 1D-array
%   containing the timestamps of the trigger events.

% Open the raw-data file, the configuration file and the timestamps file
assert(exist(filename,'file')==2,['The given file "' filename '" does not exist!']);
datafileID = fopen(filename);

config_filename = [filename '.cfg'];
assert(exist(config_filename,'file')==2,['There is no *.cfg-file ' config_filename '!']);
config = cfg2struct(config_filename);

timestamps_filename = [filename '.stamp'];
if exist(timestamps_filename,'file')==2
    timestampsfileID = fopen(timestamps_filename);
    exist_timestamps = true;
else
    disp('Warning: No timestamps file detected!');
    exist_timestamps = false;
end

% READ FILENAME TO DATA8BIT
channelnumber = config.SpectrumCard.Channel00.Enable_BOOL + config.SpectrumCard.Channel01.Enable_BOOL + config.SpectrumCard.Channel02.Enable_BOOL + config.SpectrumCard.Channel03.Enable_BOOL;
segmentsize = config.SpectrumCard.ModeSetup.Segmentsize_I32;
memsize = config.SpectrumCard.ModeSetup.Memory_I32;
number_of_recordings = memsize/segmentsize*channelnumber;
data = fread(datafileID,[segmentsize number_of_recordings], 'int8=>int8');
if channelnumber>1
    data8bit = zeros(segmentsize,memsize/segmentsize,channelnumber);
    for block=1:memsize/segmentsize
        for channel=1:channelnumber
            data8bit(:,block,channel) = data(:,(block-1)*channelnumber+channel);
        end
    end
else
    data8bit = data;
end

% READ TIMESTAMPS
% 1. When using VI dwTimestampsRead_64.vi (64bit) in LabView and big-endian ordering:
% fseek(timestampsfileID,4,'bof');
% timestamps = fread(timestampsfileID,[2*memsize/segmentsize 1],'uint64=>uint64',0,'s');
% 2. For Debugging:
% timestamps_raw = fread(timestampsfileID,[16 2*memsize/segmentsize],'uint8=>uint8');
% 3. When using VI dwTimestampsRead.vi (32bit) in LabView and big-endian ordering:
% timestamps_raw = fread(timestampsfileID,[2*memsize/segmentsize 1],'uint64=>uint64',0,'s');
% timestamps_raw(1) = mod(timestamps_raw(1),2^24); %Removing filesize bytes in the first timestamp when using 32bit timestamp reading operation in LabView
timestamps_raw = fread(timestampsfileID,[2*memsize/segmentsize 1],'uint64=>uint64');

% Removing emtpy values
timestamps = zeros(memsize/segmentsize, 1, 'uint64');
for k=1:length(timestamps)
   timestamps(k)=timestamps_raw((k-1)*2+1);
end

% CLOSE ALL OPEN FILES
fclose(datafileID);
if exist_timestamps
    fclose(timestampsfileID);
end

cd('..');
save('data8bit.mat','-v7.3','data8bit');

end

