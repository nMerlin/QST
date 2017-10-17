
for i = 56:2:90
    disp(i);
    [X1,X2,X3,piezoSign] = prepare3ChData([ num2str(i-1) '-5mW-LOonly.raw'],[num2str(i) '-5mW-LOwithDL.raw']);
    save(['2017-09-21-validated-' num2str(i) '-5mW-LOwithDL.mat'],'X1','X2','X3','piezoSign');
end