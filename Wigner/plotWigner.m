function [  ] = plotWigner( WF, varargin )
%PLOTWIGNER Plots given Wigner Function
%
% Optional Input Arguments:
%   'Image': Default is 'false'. Create an image plot of WF.
%   'Surface': Default is 'true'. Create a surface plot of WF.
%   'Narrow': Default is 'false'. Narrower axis limits.
%
% Notes:
%   p and q have to be set manually

%% Validate and parse input arguments
p = inputParser;
defaultImage = false;
addParameter(p,'Image',defaultImage,@islogical);
defaultNarrow = false;
addParameter(p,'Narrow',defaultNarrow,@islogical);
defaultSurface = true;
addParameter(p,'Surface',defaultSurface,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[image,narrow,surface] = c{:};

if narrow
    p = -6:0.125:6;
    q = p;
    nP = length(p);
    nWF = length(WF);
    shift = (nWF-nP)/2;
    WF = WF(shift+1:end-shift,shift+1:end-shift);
else
    p = -20:0.125:20;
    q = p;
end

if image
    imagesc(p,q,real(WF));
elseif surface
    h = surf(p,q,real(WF));
    set(h,'LineStyle','none');
end

end

