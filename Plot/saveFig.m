function saveFig(filename)
%SAVEFIG Wrapper for export_fig

if ~isempty(filename)
    set(gcf,'Color','w');
    [path,name,ext] = fileparts(filename);
    switch ext
        case '.jpg'
            filenameJpg = [path,name,'.jpg'];
            export_fig(sprintf('%s',filenameJpg),'-jpg','-r600');
        case '.pdf'
            filenamePdf = [path,name,'.pdf'];
            export_fig(sprintf('%s',filenamePdf),'-pdf','-painters');
        case '.fig'
            savefig(filename);
    end
end

end

