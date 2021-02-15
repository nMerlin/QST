function [ period_on, period_off, period_mean_on, period_mean_off, fliprate ] = SwitchingPhotonnumber( times, ada, g2vec, level_on, level_off, filename, varargin )
% Output: 
% period_on (s): dwell times in the "on" state
% period_off (s): dwell times in the "off" state
% period_mean_on (s): mean dwell time of the "on" state
% period_mean_off (s): mean dwell time of the "off" state
% fliprate (Hz): number of flips per time

%% Validate and parse input arguments
p = inputParser;
defaultSmooth = false; % whether data is smoothed
addParameter(p,'Smooth',defaultSmooth,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[smoothing] = c{:};
filenameOptions = ['-smooth-' num2str(smoothing) '-level_on-' num2str(level_on) '-level_off-' num2str(level_off)];

%level_on = 0.8;
%level_off = 0.2;

timesInSeconds = times;
[times,units] = convenientUnits(times,'s');

% ada_high = level_on * max(ada);
% ada_low = level_off * max(ada);
%Problem: max(ada) could be an outlier. Therefore, rather take max from
%the right peak from a histogram of ada.
% [N,edges] = histcounts(ada,'Normalisation','probability');
% adaHist = edges(1): mean(diff(edges)):edges(end-1);
% [Npks,locs] = findpeaks(N,'minPeakProminence',0.1*max(N));
% % plot the histogram of ada
% plot(adaHist, N);
% hold on;
% plot(adaHist(locs),Npks,'v');
% xlabel('Photon number');
% ylabel('Number of events');
% graphicsSettings;
% legend('Data','Found peaks','Location','best');
% savefig([filename '-histogram_photon_numbers' filenameOptions '.fig']);
% clf();
% %
% maxAda = adaHist(max(locs));
% minAda = adaHist(min(locs));
% ada_high = level_on * maxAda;
% ada_low = level_off * maxAda;

%use smoothed data to get max and min without being disturbed by outliers, and remove offset
adaSm = smooth(ada,99);
adaSm2 = adaSm - min(adaSm);
ada_high = level_on * max(adaSm2) + min(adaSm);
ada_low = level_off * max(adaSm2)  + min(adaSm);



if smoothing
    ada=adaSm;
end

t_length = length(times);
ada_digital=ada;

%Umwandlung der Daten in digitale Form
for n = 1:t_length
    if(ada(n)>ada_high)
        ada_digital(n)=1;
    elseif(ada(n)<ada_low)
        ada_digital(n)=0;
    else
        ada_digital(n)=-1;
    end
end

% plot intermediate data for diagnostic purpose
subplot(2,1,1)
plot(times,ada);
hold on;
plot(times,ada_high*ones(1,t_length),'--',times,ada_low*ones(1,t_length),'--');
ylabel('Photon number');
graphicsSettings;
h = findobj(gca,'Type','line');
set(h,'linewidth',1);
subplot(2,1,2)
plot(times, ada_digital);
xlabel(['Time (' units ')']);
graphicsSettings;
h = findobj(gca,'Type','line');
set(h,'linewidth',1);
ylim([-1.5 1.5]);
ylabel('Assignment');
savefig([filename '-intermediate-data-' filenameOptions '.fig']);
print([filename '-intermediate-data-' filenameOptions '.png'],'-dpng');
close;
%


%compute the dwell times of the "on" and "off" states
% if all photon numbers <1, there is no switching
if max(adaSm) < 1
    numberOfFlips = 0;
    period_on= 0;
    period_off = max(times);    
%     maxN-minN < range*medN
else
    period_on=zeros(1); 
    period_off=zeros(1);
    index_on = 1;
    index_off = 1;
    start_on = 0;
    start_off=0;
    numberOfFlips = 0;
    for n=1:(t_length-1)
        if( ada_digital(n+1) > ada_digital(n)&& ada_digital(n+1)==1) %switches to "on" state
            if(ada_digital(n+1)~= ada_digital(n)&& ada_digital(n)==0) %comes from "off" state
                period_end = times(n+1);
                period_off(1,index_off) = period_end - start_off;         
            end
            numberOfFlips = numberOfFlips + 1;
            index_on = index_on+1;
            start_on = times(n+1);
        elseif( ada_digital(n+1)~= ada_digital(n)&& ada_digital(n+1)==0) %switches to "off" state
            if(ada_digital(n+1)< ada_digital(n)&& ada_digital(n)==1) %comes from "on" state
                period_end = times(n+1);
                period_on(1,index_on)= period_end - start_on;            
            end
            numberOfFlips = numberOfFlips + 1;
            index_off = index_off+1;
            start_off = times(n+1);
        elseif(ada_digital(n+1)< ada_digital(n)&& ada_digital(n)==1) %comes from "on" state, goes to intermediate state
            period_end = times(n+1);
            period_on(1,index_on)= period_end - start_on;
            numberOfFlips = numberOfFlips + 1;
        elseif(ada_digital(n+1)~= ada_digital(n)&& ada_digital(n)==0) %comes from "off" state, goes to intermediate state
            period_end = times(n+1);
            period_off(1,index_off) = period_end - start_off;
            numberOfFlips = numberOfFlips + 1;
        elseif(n == t_length-1 && ada_digital(n)==1) % include the length of the last state if it is "on" until the end
            period_end = times(n);
            period_on(1,index_on)= period_end - start_on;
        elseif(n == t_length-1 && ada_digital(n)==0) % include the length of the last state if it is "off" until the end
            period_end = times(n);
            period_off(1,index_off)= period_end - start_off;

        end
    end %for n
end %if



period_mean_on = mean(period_on);
period_mean_off = mean(period_off);
fliprate = numberOfFlips / max(timesInSeconds);

%plot histograms of the dwell times
histogram(period_on,20); hold on;
histogram(period_off,20);
xlabel(['Time (' units ')']);
ylabel('Number of events');
graphicsSettings;
legend('"On" times','"Off" times','Location','best');
savefig([filename '-histogram_on_off' filenameOptions '.fig']);
print([filename '-histogram_on_off' filenameOptions '.png'],'-dpng');
f=gca; 
f.YScale = 'log';
f.XScale = 'log';
savefig([filename '-histogram_on_off_log' filenameOptions '.fig'])
print([filename '-histogram_on_off_log' filenameOptions '.png'],'-dpng');
close;

% 
% LowTimes = (times(ada < ada_low));
% LowN = (ada(ada < ada_low));
% LowG2 = (g2vec(ada < ada_low));
% 
% HighTimes = (times(ada > ada_high));
% HighN = (ada(ada > ada_high));
% HighG2 = (g2vec(ada > ada_high));

% test =[-1 1 1 -1 0 0 0 0 1 1 1 0 0 1 1 -1 -1 1 1 0 1 1 0 0 1 1 ];
% times_test = (1:26);
% index_on = 0;
% index_off = 0;
% start_off=0;
% period_on=zeros(1,2); 
% period_off=zeros(1,1);
% for n=1:(length(test)-1)
%     if( test(n+1) > test(n)&& test(n+1)==1)
%         if(test(n+1)~= test(n)&& test(n)==0)
%             period_end = times_test(n+1);
%             period_off(1,index_off) = period_end - start_off;
%         end
%         index_on = index_on+1;
%         start_on = times_test(n+1);
%     elseif( test(n+1)~= test(n)&& test(n+1)==0)
%         if(test(n+1)< test(n)&& test(n)==1)
%             period_end = times_test(n+1);
%             period_on(1,index_on)= period_end - start_on;
%         end
%         index_off = index_off+1;
%         start_off = times_test(n+1);
%     elseif(test(n+1)< test(n)&& test(n)==1)
%         period_end = times_test(n+1);
%         period_on(1,index_on)= period_end - start_on;
%     elseif(test(n+1)~= test(n)&& test(n)==0)
%         period_end = times_test(n+1);
%         period_off(1,index_off) = period_end - start_off;
%     end
% end

end




