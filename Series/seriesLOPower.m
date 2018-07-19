function seriesHDPower(varargin)
%SERIESHDPOWER creates a plot HD-Output-Power vs. LO-Power

%% Validate and parse input arguments
p = inputParser;
defaultPrepOpts = struct;
addParameter(p,'PrepOpts',defaultPrepOpts,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[prepopts] = c{:};

%% Read file names and extract LO power
filestruct = dir('raw-data/*.raw');
filestring = strjoin({filestruct.name});
fPower = regexpi(filestring,'[^ ]*mW[^ ]*.raw','match');
tokMeas = regexpi(strjoin(fPower),'\<(\d*)-','tokens');
nMeas = cellfun(@str2num,[tokMeas{:}]);
tokPow = regexpi(strjoin(fPower),'([0123456789.])*mW','tokens');
pMeas = cellfun(@str2num,[tokPow{:}]);

%% Compute signal powers
[powerHD,stdPowerHD] = deal(zeros(1,length(fPower)));
for i=1:length(fPower)
    data8bit = load8BitBinary(fPower{i},'dontsave');
    vars = var(double(data8bit));
    powerHD(i) = mean(vars);
    stdPowerHD(i) = std(vars);
end

p = polyfit(pMeas,powerHD,1);
f = polyval(p,pMeas);
errorbar(pMeas,powerHD,stdPowerHD,'o');
hold on;
plot(pMeas,f,'-');
hold off;
xlabel('LO power (mW)');
ylabel('HD Output Power (arb. unit)');
title('LO power series');
saveA5Landscape('seriesLOPower');
PowerLO = pMeas';
OutputPowerHD = powerHD';
StandardDeviation = stdPowerHD';
T = table(PowerLO,OutputPowerHD,StandardDeviation);
writetable(T,'seriesLOPower');

end

