function [data] = loadImg(filename)
% loads data from an .img file from streak camera, and saves it as an ascii
% file with .dat ending. 

fileID = fopen(filename);

%% get file information
info = fread(fileID,'int16'); % the commentary and info part is in 2 byte format.
width = info(3); % width of the image in pixels
height = info(4); % height of the image in pixels
commentlength = info(2);  %length of comment string
preLength = commentlength + 65; % length of all preliminary stuff
filetype = info(7);

%% load only the data stored in the file
frewind(fileID);
fseek(fileID,preLength-1,'bof'); % go to start of data

switch filetype
    case 0
        format = 'int8';
    case 2
        format = 'int16';
    case 3
        format = 'int32';
end
      
data = fread(fileID,[width height],format);
data = data';

fclose('all');

%% save data as ascii file
save(strrep(filename,'img','dat'),'data','-ascii');

end