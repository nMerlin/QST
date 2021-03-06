function mainWigner( nMax, mMax, varargin )
%MAINWIGNER Calculates FT(<q+1/2*q'|n>*<m|q-1/2*q'>) up to n=NMAX and m =
%MMAX. Set the discretization parameters and target directory in the source
%code.
%
% Options:
% 'parallel': Use MatLabs parallel pool for higher speed
% 'overwrite': Compute existing matrix element again

%directory = 'Z:\freefilesync-lab\matlab\QST\Wigner';
%directory = 'D:\@archive\2016-08-30-wigner-test';
directory = 'C:\Users\lab\Documents\@archived-data\Wigner';
minq=-20;
maxq=20;
qintstep=0.125; %choose 1/2^n to avoid artifacts

parallel = 0;
overwrite = 0;
if nargin > 2
    for i = 3:nargin
        eval([varargin{i-2} '=1;']);
    end
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