function selStr = selParamsToStr(selParams,varargin)
%SELPARAMSTOSTR Converts selection parameters of function 'selectRegion' to
%   human readable string for filenames.
%
% Usage:
%   selStr = selParamsToStr(selParams,varargin);
%
% Input Arguments:
%   selParams: Structure with selection parameters for 'selectRegion'.

type = selParams.Type;
pos = selParams.Position;
selStr = ['type=',selParams.Type,'-'];
switch type
    case {'fullcircle','halfcircle'}
        selStr = [selStr,'radius=',num2str(pos(1)), ...
            '-thickness=',num2str(pos(2))];
    case 'rectangle'
        selStr = [selStr,'x=',num2str(pos(1)),'-y=',num2str(pos(2)), ...
            '-width=',num2str(pos(3)),'-height=',num2str(pos(4))];
    case 'dot'
        selStr = [selStr,'x=',num2str(pos(1)),'-y=',num2str(pos(2)), ...
            '-radius=',num2str(pos(3))];
    case {'Qline','Pline'}
        selStr = [selStr,'center=',num2str(pos(1)), ...
            '-width=',num2str(pos(2))];
end

end

