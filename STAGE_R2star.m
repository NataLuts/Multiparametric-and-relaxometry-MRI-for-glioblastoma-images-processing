function [R2s, protonDensity,mask] = STAGE_R2star(mag, TE, NTE, mask)

TE = TE(:);  % a columns vector
if length(NTE)==1
    mag=mag(:,:,:,1:NTE);
    TE = TE(1:NTE);
else
    mag= mag(:,:,:,NTE);
    TE= TE(NTE);
end


dim=size(mag);
TEm=mean(TE);
TEm2=mean(TE.^2);
TE=reshape(TE,[1,1,1,dim(4)]);
TEmat=repmat(TE,[dim(1),dim(2),dim(3)]);%4D matrix

%simple linear regression
    %calculate R2s
    testzero=prod(mag,4)==0;%clean the data
    testzero=repmat(testzero,[1,1,1,dim(4)]);
    mag=log(mag);%Y
    mag(testzero==1)=0;
    clear testzero;
    magm=mean(mag,4);
    convXYr2s = mean(TEmat.*mag,4)-TEm*magm;
    varX   = TEm2 - TEm^2;
    R2s = (-1)*convXYr2s/varX;
    protonDensity  = exp(magm + R2s*TEm);
    %deltaB and phi0 maps
end