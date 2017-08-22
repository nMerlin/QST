function[O1Full, O2Full, O3Full, oThetaFull] = catMeasurements(filenameLO, startNameSIG, endNumber)
%CATMEASUREMENTS - 
%this function concatenates the orthogonal datasets from several
%measurements that have been taken after another.
%
% Input arguments: 
%   filenameLO - filename of the LO-data used for correct normalization
%   startNameSIG - filename of the raw 3-Channel-Data of the first measurement
%   endNumber - number of the last measurement
%   The names should be of the format 'number-description.raw' with
%   'description' being the same for each file.
%
% Output arguments:
% The concatenated datasets where X1 and X2 are
% orthogonal.

[startToken, description] = strtok(startNameSIG, '-');
startNumber = str2double(startToken);

rawDataContents = dir('raw-data');
name = {rawDataContents.name};
filename = cell2mat(name);

[O1Full, O2Full, O3Full, oThetaFull]  = deal([]);

for number = startNumber: endNumber
    if isempty(regexpi(filename,[num2str(number) description],'match'))
        continue;  %if the file with given number doesn't exist, continue with the next
    end        
    dispstat(['Processing Number' num2str(number)],'init','timestamp','keepthis',0);
    [X1,X2,X3, piezoSign] = prepare3ChData(filenameLO, [num2str(number) description]);
    [theta,ys] = computePhase(X1,X3,piezoSign);
    %[X1,X2,X3,theta] = segmentValidation(X1,X2,X3,ys,theta);
    [O1,O2,O3,oTheta] = selectOrthogonal(X1,X2,X3,theta,piezoSign,'Plot','hide');
    O1Full = cat(1,O1Full,O1);
    O2Full = cat(1,O2Full,O2);
    O3Full = cat(1,O3Full,O3);
    oThetaFull = cat(1,oThetaFull,oTheta);
   
end


end
