function [ HF ] = cohRegularizedPFunction(q,p,R,varargin)
%THERMHUSIMI Returns the regularized P function, with filter parameter R,  for a coherent state
%   Detailed explanation goes here

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultP0 = 1; 
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,p0,q0] = c{:};

p = p';
distance = sqrt(((q-q0).^2 + (p-p0).^2)/(2*norm)^2);
HF = (besselj(1,2*R*distance)./(sqrt(pi)*distance)).^2;
HF = HF./sum(sum(HF));

% The Glauber-Sudarshan P-distribution of the coherent state is a delta
% distribution delta(alpha - alpha_0). Thus, the regularized P function of the coherent state is the
% filter function Omega(alpha - alpha_0).


end

