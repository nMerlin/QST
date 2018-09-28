function [ HF ] = thermHusimi(q,p,nPhotons,varargin)
%THERMHUSIMI Returns Q(q,p) for a thermal state with NPHOTONS.
%   Detailed explanation goes here

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm] = c{:};

disc = mean(diff(q));
HF = zeros(length(q),length(p));
for iP = 1:length(p)
    HF(:,iP)=disc^2*(1/norm)^2*1/(pi*(nPhotons+1))* ...
        exp(-((q*norm).^2+(p(iP)*norm).^2)/(nPhotons+1));
end
%HF = HF./sum(sum(HF)); % Usually not necessary

end

