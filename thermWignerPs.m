function WF = thermWignerPs(q,p,qps,pps,r,nPhotons,varargin)
%THERMWIGNERPS Returns the conditional Wigner function W(q,p,qps,pps) for
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
WF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    WF(:,iP) = disc^2*norm^2*2*(1+nPhotons*t^2)/ ...
        (pi*(1+nPhotons+nPhotons*r^2))*exp(-2*((1+nPhotons*t^2)*norm* ...
        qrot+nPhotons*r*t*qps).^2/((1+nPhotons*t^2)* ...
        (1+nPhotons+nPhotons*r^2))).*exp(-2*((1+nPhotons*t^2)*norm* ...
        prot+nPhotons*r*t*pps).^2/((1+nPhotons*t^2)* ...
        (1+nPhotons+nPhotons*r^2)));
end
%WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end