function [ data8bit, config, timestamps ] = load8BitBinary( filename )
%LOAD8BITBINARY Loads 8bit binary datafiles, the configuration file and the
%timestamp file for a single multiple recording measurement.
%
%   DATA8BIT = LOAD8BITBINARY(filename) Returns a 2D-Array...

% Open the raw-data file, the configuration file and the timestamps file
assert(exist(filename,'file')==2,['The given file "' filename '" does not exist!']);
datafileID = fopen(filename);
data8bit = [];

config_filename = [filename '.cfg'];
assert(exist(config_filename,'file')==2,['There is no *.cfg-file ' config_filename '!']);
config = cfg2struct(config_filename);

timestamps_filename = [filename '.stamp'];
timestamps = [];
if exist(timestamps_filename,'file')==2
    timestampsfileID = fopen(timestamps_filename);
    exist_timestamps = true;
else
    disp('Warning: No timestamps file detected!');
    exist_timestamps = false;
end

% Determine number of channels
% config.SpectrumCard.Channel01.Enable

% Close all open files
fclose(datafileID);
if exist_timestamps
    fclose(timestampsfileID);
end

end

