function [  ] = plotCorrelation( j, data, window )
%PLOTCORRELATION
X = 1:floor(window);
Y = zeros(floor(window),1);
for detuning=1:floor(window)
    Y(detuning) = correlation(j, data, detuning, window);
end

plot(X,Y);

end

