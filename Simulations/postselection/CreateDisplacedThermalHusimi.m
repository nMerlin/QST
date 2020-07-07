function [ HusFunc,xx,yy ] = CreateDisplacedThermalHusimi( n, Q0,P0,WLength, WRes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


QuadVals=-WLength:WRes:WLength;
MatSize=round(2*(WLength/WRes)+1);

HusFunc=zeros(MatSize,MatSize);
[xx,yy] = meshgrid(QuadVals,QuadVals);

HusFunc=CreateDisplacedThermalHusimiAtQuad( n,Q0,P0,xx,yy,WRes );

%This defintion of the coherent state means that the quads do not correspond to
%the amplitudes, <Qx^2>=<Qy^2>=0.5 for the vacuum state

%Should be 1, but just to make sure:
HusFunc=HusFunc./sum(sum(HusFunc));

end

