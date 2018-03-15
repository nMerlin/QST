function plotG2(times,values)
%PLOTG2 Plot g2 values against recording times (expects SI units)
%
% Input Arguments:
%   times: vector with time axis
%   values: vector with g2 values

%% Create plot
[times,units] = convenientUnits(times,'s');
plot(times,values);
axis tight;
ylim([0,3]);
xlabel(['Time (',units,')']);
ylabel('g2(0)');

end

