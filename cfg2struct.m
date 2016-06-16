function Struct = ini2struct(FileName)
% Parses .ini file
% Returns a structure with section names, subsection names and keys as fields.
% 
% Based on:
% freeb (2014). ini2struct
% (http://www.mathworks.com/matlabcentral/fileexchange/45725-ini2struct/content/ini2struct.m),
% MATLAB Central File Exchange. Retrieved June 15, 2016.
% 
% Modified by:
% Johannes Thewes
% June 15, 2016

f = fopen(FileName,'r');                    % open file
while ~feof(f)                              % and read until it ends
    s = strtrim(fgetl(f));                  % remove leading/trailing spaces
    if isempty(s) || s(1)==';' || s(1)=='#' % skip empty & comments lines
        continue
    end
    if s(1)=='['                            % section header
        if length(strfind(s,'.'))==1        % if there is section and subsection
            subsection = true;
            C = strsplit(s,'.');
            s1 = C{1};
            s2 = C{2};
            Section = genvarname(strtok(s1(2:end),']'))
            Subsection = genvarname(strtok(s2(1:end),']'));
            Struct.(Section).(Subsection) = [];
        else
            subsection = false;
            Section = genvarname(strtok(s(2:end), ']'));
            Struct.(Section) = [];              % create field
        end
        continue
    end
    
    [Key,Val] = strtok(s, '=');             % Key = Value ; comment
    Val = strtrim(Val(2:end));              % remove spaces after =
    
    if isempty(Val) || Val(1)==';' || Val(1)=='#' % empty entry
        Val = [];
    elseif Val(1)=='"'                      % double-quoted string
        Val = strtok(Val, '"');
    elseif Val(1)==''''                     % single-quoted string
        Val = strtok(Val, '''');
    else
        Val = strtok(Val, ';');             % remove inline comment
        Val = strtok(Val, '#');             % remove inline comment
        Val = strtrim(Val);                 % remove spaces before comment
        
        [val, status] = str2num(Val);       %#ok<ST2NM>
        if status, Val = val; end           % convert string to number(s)
    end
    
    if ~exist('Section', 'var')             % No section found before
        Struct.(genvarname(Key)) = Val;
    elseif subsection
        Struct.(Section).(Subsection).(genvarname(Key)) = Val;
    else                                   % Section found before, fill it
        Struct.(Section).(genvarname(Key)) = Val;
    end

end
fclose(f);
