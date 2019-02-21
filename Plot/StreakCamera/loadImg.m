function [data, time] = loadImg(filename)
% loads data from an .img file from streak camera, and saves it as an ascii
% file with .dat ending. 
% time contains the time scale in ps (switch case for unit?). 

fileID = fopen(filename);

%% get file information
info = fread(fileID,'int16'); % the commentary and info part is in 2 byte format.
width = info(3); % width of the image in pixels
height = info(4); % height of the image in pixels
commentlength = info(2);  %length of comment string
preLength = commentlength + 65; % length of all preliminary stuff
filetype = info(7);

%get comment string and start of time information
frewind(fileID);
fseek(fileID,65,'bof');
comment = fread(fileID,commentlength,'*char')';
starttimeToken = regexpi(comment,'ScalingYScalingFile="#([0123456789]*),','tokens');
startTime = str2double(cell2mat(starttimeToken{1}));

%% load only the data stored in the file
frewind(fileID);
fseek(fileID,preLength-1,'bof'); % go to start of data

switch filetype
    case 0
        fmt = 'int8';
    case 2
        fmt = 'int16';
    case 3
        fmt = 'int32';
end
      
data = fread(fileID,[width height],fmt);
data = data';

% get the time scale
frewind(fileID);
fseek(fileID,startTime,'bof');
time = fread(fileID,height,'single');

fclose('all');

%% save data as ascii file
save(strrep(filename,'img','dat'),'data','-ascii');
save(strrep(filename,'.img','-time.txt'),'time','-ascii');

end