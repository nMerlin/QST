function plotPulsePostselection(X, channel)
%This function makes a plot showing the dependence between pulses shifted 
%about delay. Each pulse is normalized so the bins are comparable.
%The values of one pulse are sorted into bins and than a pulse shifted
%about delay is posteselected accordingly and the mean values computed. 
%If the two pulses are independent, they should be zero. 
%This is done for many pulses shifted about a certain delay and averaged.

varBins = 10;
maxDelay = 5;

%% Preparation

XCh = X(:,:,channel);

for i = 1:size(XCh,1)
XCh(i,:) = XCh(i,:) / max(abs(XCh(i,:))); %Normalize each pulse
end

%% Sorting of Pulse Quadrature Values into Bins 
meanXBinned = NaN(size(XCh,1)-maxDelay,varBins);
for delay = 1:maxDelay
    for pulse = 1:size(XCh,1)-delay
    % Sorting of Pulse Quadrature Values into Bins 
        [N,varEdges,bin] = histcounts(XCh(pulse,:),varBins);
        [~,I] = sort(bin);
      % Sorting the values of the pulse shifted about delay into these bins  
        Xdelay = XCh(pulse+delay-1,I);
        XOut = NaN(max(N), varBins);
        for iInterval = 1 : varBins
            start = 1+sum(N(1:iInterval-1));
            stop = start+N(iInterval)-1;
            XOut(1:N(iInterval),iInterval) = Xdelay(start:stop);
        end
        meanXBinned(pulse,:) = mean(XOut, 'omitnan');
    end   
    meanMeanXBinned = mean(meanXBinned, 'omitnan'); %Average over all pulses with this delay
    xAxis = varEdges(1:end-1)+min(diff(varEdges))/2;
    plot(xAxis,meanMeanXBinned,'-o','linewidth',2, 'color',[1 0 0]*delay/maxDelay,...
        'DisplayName',['Delay ' num2str(delay-1)]);
    hold on;
end

legend('location','best');
xlabel('Quadrature Bins');
ylabel('Mean of Postselected Quadratures');
title(['Binned Quadrature Values (',num2str(varBins),' Bins)']);
axis([-1.5 1.5 -1.5 1.5]);

end