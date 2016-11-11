function [  ] = plotWigner( WF, varargin )
%PLOTWIGNER Creates a plot of the given Wigner Function
%
%   Noes:
%   p and q have to be set manually

p = -20:0.125:20;
q = -20:0.125:20;

image = 0;
surface = 0;
narrow = 0;
if nargin > 1
    for i = 2:nargin
        eval([varargin{i-1} '=1;']);
    end
end

if narrow == 1
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

if image == 1
    imagesc(p,q,real(WF));
elseif surface == 1
    h = surf(p,q,real(WF));
    set(h,'LineStyle','none');
end

end

