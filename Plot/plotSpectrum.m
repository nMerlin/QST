function plotSpectrum(data,varargin)
%PLOTSPECTRUM plots an optical spectrum
%
%   Input Arguments:
%       data: array with two columns; first column contains x-values and
%           second column contains y-values

plot(data(:,1),data(:,2),'LineWidth',2);
axis tight
set(gca,'YLim',[0 65000]);
ylabel('Counts');
xlabel('Wavelength (nm)');
title('Spectrum');

end

