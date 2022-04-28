function [ HusFunc ] = CreateCoherentHusimi( x0,y0,HLength,HRes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


QuadVals=-HLength:HRes:HLength;
MatSize=2*(HLength/HRes)+1;

HusFunc=zeros(MatSize,MatSize);
[xx,yy] = meshgrid(QuadVals,QuadVals);

HusFunc=CreateCoherentHusimiAtQuad( x0,y0,xx,yy,HRes );



%This defintion of the coherent state means that
% <Qx^2>=<Qy^2>=1
HusFunc=HusFunc./sum(sum(HusFunc));

end

