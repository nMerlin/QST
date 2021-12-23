function [] = doubleslitSimulation()
%Distance between slits in meter
dSlit = 133e-6;
%Wavelength in m 
lambda = 832e-9;
% length of free space path in meter
D = 0.6;

%Intensity ratio between the slits:
ratio = 0:0.1:1;
ratioMatrix = ones(1,1, length(ratio));
ratioMatrix(1,1,:) =  ratio;
%phase difference between the two light sources:
DeltaTheta = 0:0.1:pi/2;
DeltaThetaMatrix = ones(1,1,1, length(DeltaTheta));
DeltaThetaMatrix(1,1,1,:) =  DeltaTheta;

%% 
%original phases at the beginning:
theta1 = 0:0.01:2*pi;  
theta1 = theta1';
theta2 = theta1 + DeltaThetaMatrix;

%horizontal Position after the free space path (zero is middle) in meter:
position = -0.01:0.0001:0.01;

%Path lengths from doubleslit to measurement points for slit 1 (left from
%the middle):
s1 = sqrt((position + dSlit/2).^2 + D^2);
s2 = sqrt((position - dSlit/2).^2 + D^2); 

%Phases from the path:
phase1 = mod(2*pi/lambda * s1, 2*pi);
phase2 = mod(2*pi/lambda * s2, 2*pi);

% Intensity pattern:

Int = abs( sqrt(ratioMatrix).*cos(theta1 + phase1)  +  cos(theta2 + phase2) ).^2;

%averaging over all phases of one light source 
%(corresponding to averaging over one cycle of the electrical field (with fixed phase difference
%between the light sources):

Int = mean(Int,1);

%plot all intensity ratios for one DeltaTheta:
IntReshaped = Int(:,:,:,1);
IntReshaped = reshape(IntReshaped,[size(IntReshaped,2) size(IntReshaped,3)]);
plot(IntReshaped);
plot(position*1000,IntReshaped, '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '-DeltaTheta-' num2str(DeltaTheta(1)) '.fig']);
clf();

%averaging over the intensity ratios:
IntAveraged = mean(Int,3);
plot(position*1000,IntAveraged(1,:,1,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '-DeltaTheta-' num2str(DeltaTheta(1)) '.fig']);
clf();

%plot all deltaThetas for intensity ratio 1:
IntReshaped = Int(:,:,end,:);
IntReshaped = reshape(IntReshaped,[size(IntReshaped,2) size(IntReshaped,4)]);
plot(IntReshaped);
plot(position*1000,IntReshaped, '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-IntensityRatio-' num2str(ratio(end)) '.fig']);
clf();

%averaging over the deltaThetas for one intensity ratio:
IntAveraged = mean(Int,4);
plot(position*1000,IntAveraged(1,:,end,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-IntensityRatio-' num2str(ratio(end)) '.fig']);
clf();

%averaging over deltaThetas and intensity ratios:
IntAveraged = mean(Int,4);
IntAveraged = mean(IntAveraged,3);
plot(position*1000,IntAveraged(1,:,1,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-AveragedIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '.fig']);
clf();


end