function HF = thermHusimiPs(q,p,qps,pps,r,nPhotons,varargin)
%THERMHUSIMIPS Returns the conditional Husimi-Q function H(q,p,qps,pps) for
%a thermal state with NPHOTONS photons. r is the reflection coefficient of
%the first signal beamsplitter in a 12-port homodyne detector into the
%target arm.
%
% Formulation is in alpha-space but choice of normalization has to be taken
% into account.

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultTheta = 0; % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,theta] = c{:};

t = sqrt(1-r^2);
disc = mean(diff(q));
HF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    HF(:,iP) = disc^2*norm^2*(1+nPhotons*t^2)/(pi*(1+nPhotons))* ...
        exp((qps^2+pps^2)/(1+nPhotons*t^2)).* ...
        exp(-(r*qps+t*norm*qrot).^2-(r*pps+t*norm*prot).^2).* ...
        exp(-1/(1+nPhotons)*((t*qps-r*norm*qrot).^2+ ...
        (t*pps-r*norm*prot).^2));
end
%WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end
