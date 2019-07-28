% Function handle: Rotation
% Purpose: This function is used to rotate the parts of the gyroscope
% Inputs: 
% xi: Initial x-coordinates
% yi: Initial y-coordinates
% zi: Initial z-coordinates
% R: Rotation matrix chosen based on which frame the object is in
% Outputs: 
% xf: Rotated x-coordinates 
% yf: Rotated y-coordinates
% zf: Rotated z-coordinates

function [xf,yf,zf]=rotate(xi,yi,zi,R)

I=size(xi,1);
J=size(xi,2);

xf=zeros(I,J);
yf=zeros(I,J);
zf=zeros(I,J);

for ii=1:I
    for jj=1:J
        vector=[xi(ii,jj);yi(ii,jj);zi(ii,jj)];
        vector=R*vector;
            xf(ii,jj)=vector(1);
            yf(ii,jj)=vector(2);
            zf(ii,jj)=vector(3);
    end
end