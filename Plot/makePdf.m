function makePdf(figFile,pdfpath)
%MAKEPDF Exports the *.fig file in 'figFile' to a pdf with the same name in
%the 'pdfpath' directory using 'export_fig'.

fig = openfig(figFile);
[~,name] = fileparts(figFile);
filenamePdf = [pdfpath,name,'.pdf'];
export_fig(sprintf('%s',filenamePdf),'-pdf','-painters');
close(fig);

end

