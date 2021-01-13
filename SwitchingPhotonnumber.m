function [ period_on, period_off ] = SwitchingPhotonnumber( filename, level_on,level_off )



load(filename, 'times', 'ada', 'g2vec');

%ada=smooth(ada,99);
t_length = length(times);
ada_high = level_on * max(ada);
ada_low = level_off * max(ada);
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


%Berechnung der Zeitspannen der on und off Zustände
period_on=zeros(1,10); 
period_off=zeros(1,10);
index_on = 0;
index_off = 0;
start_off=0;
for n=1:(t_length-1)
if( ada_digital(n+1) > ada_digital(n)&& ada_digital(n+1)==1)
    if(ada_digital(n+1)~= ada_digital(n)&& ada_digital(n)==0)
        period_end = times(n+1);
        period_off(1,index_off) = period_end - start_off;
    end
    index_on = index_on+1;
    start_on = times(n+1);
elseif( ada_digital(n+1)~= ada_digital(n)&& ada_digital(n+1)==0)
    if(ada_digital(n+1)< ada_digital(n)&& ada_digital(n)==1)
        period_end = times(n+1);
        period_on(1,index_on)= period_end - start_on;
    end
    index_off = index_off+1;
    start_off = times(n+1);
elseif(ada_digital(n+1)< ada_digital(n)&& ada_digital(n)==1)
    period_end = times(n+1);
    period_on(1,index_on)= period_end - start_on;
elseif(ada_digital(n+1)~= ada_digital(n)&& ada_digital(n)==0)
    period_end = times(n+1);
    period_off(1,index_off) = period_end - start_off;
end
end



period_mean_on = mean(period_on);
period_mean_off = mean(period_off);


%Darstellung der Perioden in Histogram
histogram(period_on, 100)
xlabel('time(s)');
savefig('histogram_on.fig')
f=gca; 
f.YScale = 'log';
f.XScale = 'log';
savefig('histogram_on_log.fig')
print('histogram_on.png','-dpng');

histogram(period_off, 100)
xlabel('time(s)');
f=gca; 
f.YScale = 'log';
f.XScale = 'log';
savefig('histogram_off_log.fig')
print('histogram_off.png','-dpng');



clear;
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




