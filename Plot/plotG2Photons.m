function plotG2Photons(times,values)
%PLOTG2PHOTONS Plot photon numbers of a g2 measurement against time axis
%
% Input Arguments:
%   times: vector with time axis (in seconds)
%   values: vector with photon number values

%% Create plot
[times,units] = convenientUnits(times,'s');
plot(times,values);
axis tight;
ylim([0,300]);
xlabel(['Time (',units,')']);
ylabel('# Photons');

end

