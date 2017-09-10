function plotHusimi(O1,O2,iSelect)
%PLOTHUSIMI Plot histogram of Husimi-Q function
    [H, binsO1, binsO2] = histogram2D(O1,O2);
    imagesc(binsO1,binsO2,H); axis on; colormap hot; hold on;
    plot(O1(iSelect),O2(iSelect),'.'); hold off;
    set(gca,'XLim',[-7,7],'YLim',[-7,7]);
    pbaspect([1 1 1]); % ensure that a circle doesn't look elliptic
    title('Histogram of Husimi-Q and Selected Region');
    xlabel('X1');
    ylabel('X2');
end

