function plotHusimi(O1,O2,iSelect)
%PLOTHUSIMI Plot histogram of Husimi-Q function
    [H, binsO1, binsO2] = histogram2D(O1,O2);
    imagesc(binsO1,binsO2,H); axis on; colormap hot; hold on;
    plot(O1(iSelect),O2(iSelect),'.'); hold off;
end

