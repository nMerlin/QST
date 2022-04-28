function [] = concatenateG2Data(varargin)
%this function concatenates g2data that were measured at the
%same power and saves them in one file. It also plots them in time domain
%and fft. 

%% Validate and parse input arguments
p = inputParser;
defaultX = 'X1'; % which quadrature channel should be used
addParameter(p,'UseX',defaultX);
defaultFolder = 'g2data-'; % which folder contains the data
addParameter(p,'Folder',defaultFolder);
defaultPowersteps = 150:20:350; % the powers in mW that were used 
addParameter(p,'Powersteps',defaultPowersteps);
parse(p,varargin{:});
c = struct2cell(p.Results);
[folder,powersteps,useX] = c{:};

folder = [folder useX];
[filenames,numbers,Is]= getParametersFromFilenames('Folder',folder);
stepwidth = mean(diff(powersteps));
if ~exist(['g2data-' useX '-Concatenated'],'dir')
    mkdir(['g2data-' useX '-Concatenated'])
end

for powerstep = powersteps
  % powers = Is((Is>(powerstep - stepwidth/2)) & (Is< (powerstep + stepwidth/2)) );
   files = filenames((Is>(powerstep - stepwidth/2)) & (Is< (powerstep + stepwidth/2)) ); 
   adaCon = [];
   g2Con = [];
   timesCon = [];
   for fileI = 1:length(files)
      load([folder '\' cell2mat(files(fileI))], 'times', 'ada', 'g2vec'); 
      adaCon = [adaCon ada];
      g2Con = [g2Con g2vec'];
      if fileI >1
          times = times + timesCon(end);
      end
      timesCon = [timesCon times];
   end 
   ada = adaCon;
   g2vec = g2Con;
   times = timesCon;
   save(['g2data-' useX '-Concatenated' '\' cell2mat(files(1)) '-Concatenated.mat'],'ada','g2vec','times');
   plotG2vsTime(timesCon, g2Con, adaCon, ['Plots-' useX '\' strrep(cell2mat(files(1)),'.mat','-Concatenated')]);
   plotFFTfromNphotonsVsTimes(adaCon,timesCon,['Plots-' useX '\' strrep(cell2mat(files(1)),'.mat','-Concatenated')],...
       'Axis',[0 1000 0 1]);
   clf();
end
end 