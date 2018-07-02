function WF = fockWigner(q,p,nPhotons,varargin)
%FOCKWIGNER Returns W(q,p) for a fock state with NPHOTONS photons
%
% Formula from:
% https://www.qms.physik.uni-rostock.de/fileadmin/uni-rostock/Alle_MNF/Physik_Qms/Lehre_Scheel/quantenoptik/Quantenoptik-Vorlesung4.pdf

%% Validate and parse input arguments
parser = inputParser;
defaultSigma = 1/sqrt(2); % Select another norm if applicable
addParameter(parser,'Sigma',defaultSigma,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[sigma] = c{:};

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    WF(:,iP) = 1/(2*pi*sigma^2)*(-1)^nPhotons*exp(-1/(2*sigma^2)* ...
        (q.^2+p(iP).^2)).*laguerreL(nPhotons,1/sigma^2*(q.^2+p(iP).^2));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end