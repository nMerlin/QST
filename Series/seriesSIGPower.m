function seriesSIGPower(varargin)
%SERIESSIGPOWER creates a plot HD-Output-Power vs. SIG-Power

%% Validate and parse input arguments
p = inputParser;
defaultConfigFile = '';
addParameter(p,'ConfigFile',defaultConfigFile,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[configfile] = c{:};

%% Read file names and extract power figure
filestruct = dir('raw-data/*.raw');
filestring = strjoin({filestruct.name});
fPower = regexpi(filestring,'[^ ]*mW[^ ]*.raw','match');
tokMeas = regexpi(strjoin(fPower),'\<(\d*)-','tokens');
nMeas = cellfun(@str2num,[tokMeas{:}]);
tokPow = regexpi(strjoin(fPower),'([0123456789.])*mW','tokens');
pMeas = cellfun(@str2num,[tokPow{:}]);

%% Load config file
if ~strcmp(configfile,'')
    conf = csvread(configfile);
    pMeas = conf(:,2);
end

%% Compute signal powers
[powerHD,stdPowerHD] = deal([]);
for i=1:length(fPower)
    if ~strcmp(configfile,'') && ~ismember(i,conf(:,1))
        continue;
    end
    [data8bit,~,~] = load8BitBinary(fPower{i},'dontsave');
    vars = var(double(data8bit));
    powerHD(end+1,1) = mean(vars);
    stdPowerHD(end+1,1) = std(vars);
end

p = polyfit(pMeas,powerHD,1);
f = polyval(p,pMeas);
errorbar(pMeas,powerHD,stdPowerHD,'o');
hold on;
plot(pMeas,f,'-');
hold off;
xlabel('SIG power (mW)');
ylabel('HD Output Power (arb. unit)');
title('SIG power series');
[~,configname,~] = fileparts(configfile);
saveA5Landscape(['seriesSIGPower',configname]);
PowerSIG = pMeas;
OutputPowerHD = powerHD;
StandardDeviation = stdPowerHD;
T = table(PowerSIG,OutputPowerHD,StandardDeviation);
writetable(T,['seriesSIGPower',configname]);

end

