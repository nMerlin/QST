function mainWigner( nMax, mMax, varargin )
%MAINWIGNER Calculates FT(<q+1/2*q'|n>*<m|q-1/2*q'>) up to n=NMAX and m =
%MMAX. Set the discretization parameters and target directory in the source
%code.
%
% Options:
% 'parallel': Use MatLabs parallel pool for higher speed
% 'overwrite': Compute existing matrix element again
%resulting matrix will be saved in the given DIRECTORY. The function will
%be discretized in phase space from p,q = MINQ to MAXQ in steps of
%QINTSTEP.
%% Validate and parse input arguments
p = inputParser;
defaultDirectory = 'C:\Users\lab\Documents\@archived-data\Wigner';
%directory = 'Z:\freefilesync-lab\matlab\QST\Wigner';
%directory = 'D:\@archive\2016-08-30-wigner-test';
%directory = 'C:\Users\lab\Documents\@archived-data\Wigner-Resolution-0.25';
addParameter(p,'Directory',defaultDirectory,@isstr);
defaultMinq = -20;
addParameter(p,'Minq',defaultMinq,@isnumeric);
defaultMaxq = 20;
addParameter(p,'Maxq',defaultMaxq,@isnumeric);
defaultQintstep = 0.125;  %choose 1/2^n to avoid artifacts
addParameter(p,'Qintstep',defaultQintstep,@isnumeric);
defaultParallel = 0; %choose 1 for parallel processing
addParameter(p,'Parallel',defaultParallel,@isnumeric);
defaultOverwrite = false;
addParameter(p,'Overwrite',defaultOverwrite,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[directory,maxq,minq,overwrite,parallel,qintstep] = c{:};

%% Create directory
if ~exist([pwd directory],'dir')
    mkdir(directory)
end

dispstat('','init');
dispstat('Beginning calculation...','keepthis','timestamp');

if parallel==1
    parfor i=0:(nMax+1)*(mMax+1)-1
        nVal = floor(i/(mMax+1));
        mVal = mod(i,mMax+1);
        calcWignerTable(nVal, mVal, directory, minq, maxq, qintstep, ...
            parallel, overwrite);
        dispstat('','init');
        dispstat(strcat('FT(<q+1/2*q''|',int2str(nVal),'>*<', ...
            int2str(mVal),'|q-1/2*q''>)',' computed!'), ...
            'timestamp','keepthis');
    end
else
    for i=0:(nMax+1)*(mMax+1)-1
        nVal = floor(i/(mMax+1));
        mVal = mod(i,mMax+1);
        dispstat(strcat('Computing',' FT(<q+1/2*q''|', ...
            int2str(nVal),'>*<',int2str(mVal),'|q-1/2*q''>)','...'), ...
            'timestamp','keepthis');
        calcWignerTable(nVal, mVal, directory, minq, maxq, qintstep, ...
            parallel, overwrite);
    end
end

dispstat('Calculation finished!','timestamp');

end