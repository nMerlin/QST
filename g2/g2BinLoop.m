function [] = g2BinLoop(X, BinVector, ResFix, ResVector, BinFix)
%This function evaluates the function g2Bins for various amounts of bins
%and various resolutions. 
%
% Input Arguments:
%   X - Quadratures from a continuous quantum state measurement (equal time
%       spacing between all points)
%   BinVector: Contains the amounts of Bins to be used, while using a fixed
%   resolution ResFix
%   ResVector: Contains the resolutions to be used, while using a fixed
%   amounnt of bins BinFix


% Make a resolution vector which is logarithmic over all frequencys:
% f = 0:0.3:17;
% g = exp(f);
% ResVector = 75.4e6 ./ g;
% ResVector = round(ResVector);

%% loop over different amounts of Bins 
meang2VectorBins = zeros(size(BinVector));
 i = 1;
for varBins = BinVector 
    [~, ~, meang2] = g2Bins(X, ResFix, varBins);
    meang2VectorBins(i) = meang2;
    i = i + 1;
end

%% loop over different initial resolutions 
meang2VectorRes = zeros(size(ResVector));
 i = 1;
for Resolution = ResVector 
    [~, ~, meang2] = g2Bins(X, Resolution, BinFix);
    meang2VectorRes(i) = meang2;
    i = i + 1;
end

fVector = 75.4e6 ./ ResVector;

%% Plot
subplot(3,1,1);
plot(BinVector,meang2VectorBins,'o');
%ylim([-0.5 3]);
xlabel('Number of Bins');
ylabel('mean g2 in the middle of range');
text(0.5, 0.75, ['Fixed Resolution: ' num2str(ResFix,'%.0f')],'Units','normalized');
hold on;
subplot(3,1,2);
%plot(ResVector,meang2VectorRes,'o');
semilogx(ResVector,meang2VectorRes,'o');
xlabel('Initial Resolution');
ylabel('mean g2 in the middle of range');
text(0.5, 0.75, ['Fixed Number of Bins: ' num2str(BinFix,'%.0f')],'Units','normalized');
subplot(3,1,3);
%plot(ResVector,meang2VectorRes,'o');
semilogx(fVector,meang2VectorRes,'o');
xlabel('Frequency (Hz)');
ylabel('mean g2 in the middle of range');
text(0.5, 0.75, ['Fixed Number of Bins: ' num2str(BinFix,'%.0f')],'Units','normalized');
savefig('g2-Binned-loop.fig');
save('BinSeriesData.mat','BinVector','meang2VectorBins','ResVector','ResVector','meang2VectorRes');
clf();

end