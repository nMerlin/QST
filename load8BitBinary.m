function [ data8bit, config, timestamps ] = load8BitBinary( filename )
%LOAD8BITBINARY Loads 8bit binary datafiles, the configuration file and the
%timestamp file for a single multiple recording measurement.
%
%   DATA8BIT = LOAD8BITBINARY(filename) Returns a 2D-Array...

assert(exist(filename,'file')==2,['The given file "' filename '" does not exist!']);
datafileID = fopen(filename);

config_filename = [filename '.cfg'];
assert(exist(config_filename,'file')==2,['There is no *.cfg-file ' config_filename '!']);
configfileID = fopen(config_filename);

fclose(datafileID);
fclose(configfileID);
end

