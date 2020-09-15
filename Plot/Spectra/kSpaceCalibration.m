function [] = kSpaceCalibration(filename)
%% load data
    %cd('raw-data');
    data = textread(filename,'','delimiter',',','headerlines',1);
    %data = textread(filename,'','headerlines',1);
    X = data(:,1); % wavelength
    Y = data(:,2); % pixel position
    Int = data(:,3); %intensity
    
 
 %% sort data
if Y(1) == 0
    Y=Y+1;
end
 n = length(Y)/max(Y);
 X = X(1:n);
 Y = 1:max(Y);
 Y = Y';
 Int = reshape(Int, [n max(Y)]);
 Int = Int';
 
surf(X, Y, log(Int));
colorbar;
view(180,-90);
shading flat;
axis tight;
 axis equal;


savefig(strrep(filename,'csv','fig'));
end