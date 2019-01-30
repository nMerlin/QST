function [ PF ] = thermGlauber(q,p,nPhotons,varargin)
%THERMHUSIMI Returns P(q,p) for a thermal state with NPHOTONS.

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm] = c{:};

disc = mean(diff(q));
PF = zeros(length(q),length(p));
for iP = 1:length(p)
    PF(:,iP)=disc^2*norm^2*1/(pi*nPhotons)* ...
        exp(-((q*norm).^2+(p(iP)*norm).^2)/nPhotons);
end
%PF = PF./sum(sum(HF)); % Usually not necessary

end