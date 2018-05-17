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
powerSIG = [];
for i=1:length(fPower)
    if ~strcmp(configfile,'') && ~ismember(i,conf(:,1))
        continue;
    end
    [data8bit,~,~] = load8BitBinary(fPower{i},'dontsave');
    powerSIG(end+1,1) = mean(var(double(data8bit)));
end

p = polyfit(pMeas,powerSIG,1);
f = polyval(p,pMeas);
plot(pMeas,powerSIG,'o',pMeas,f,'-');
xlabel('SIG power (mW)');
ylabel('HD Output Power (arb. unit)');
title('SIG power series');
[~,configname,~] = fileparts(configfile);
saveA5Landscape(['seriesSIGPower',configname]);
PowerLO = pMeas';
OutputPowerSIG = powerSIG';
T = table(PowerLO,OutputPowerSIG);
writetable(T,['seriesSIGPower',configname]);

end

