function selParams = selStrToParams(selStr,varargin)
%selStrToParams Conerts human readable string with postselection parameters
%used in 'selectRegion' to structure with these parameters.

type = regexp(selStr,'type=(\w*)','tokens');
if ~isempty(type)
    type = type{1}{1};
else
    type = '';
end
switch type
    case {'fullcircle','halfcircle'}
        r = regexp(selStr,'radius=([0-9.]*)','tokens');
        pos(1) = str2double(r{1}{1});
        t = regexp(selStr,'thickness=([0-9]*.[0-9]*)','tokens');
        pos(2) = str2double(t{1}{1});
    case 'rectangle'
        x = regexp(selStr,'x=([0-9.]*)','tokens');
        pos(1) = str2double(x{1}{1});
        y = regexp(selStr,'y=([0-9.]*)','tokens');
        pos(2) = str2double(y{1}{1});
        w = regexp(selStr,'width=([0-9.]*)','tokens');
        pos(3) = str2double(w{1}{1});
        h = regexp(selStr,'height=([0-9]*.[0-9]*)','tokens');
        pos(4) = str2double(h{1}{1});
    case 'dot'
        x = regexp(selStr,'x=([0-9.]*)','tokens');
        pos(1) = str2double(x{1}{1});
        y = regexp(selStr,'y=([0-9.]*)','tokens');
        pos(2) = str2double(y{1}{1});
        r = regexp(selStr,'radius=([0-9]*.[0-9]*)','tokens');
        pos(3) = str2double(r{1}{1});
    case {'Qline','Pline'}
        c = regexp(selStr,'center=([0-9.]*)','tokens');
        pos(1) = str2double(c{1}{1});
        w = regexp(selStr,'width=([0-9]*.[0-9]*)','tokens');
        pos(2) = str2double(w{1}{1});
end
if ~isempty(type)
    selParams.Type = type;
    selParams.Position = pos;
else
    selParams = struct('Type',{},'Position',{});
end

end

