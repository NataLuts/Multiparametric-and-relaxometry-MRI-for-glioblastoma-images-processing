%by Saifeng Liu
%June 20,2012
function [hpima,mask]=hpfilter_swim(phase,sizex,sizey, dx, dy)
%homodyne high-pass filter
dim=size(phase);
%change the size of the hpfilter according to the dimension of the image
sizex=sizex*dim(1)*dx/512/0.5;
sizey=sizey*dim(2)*dy/512/0.5;
sizex=round(sizex);
sizey=round(sizey);
if mod(sizex,2)==1
    sizex=sizex+1;
end
if mod(sizey,2)==1
    sizey=sizey+1;
end
mag=ones(dim);
hpima=zeros(dim);
t1=1:sizex;
t2=1:sizey;
f1=0.5*(1+cos(2*pi*(t1-sizex/2-1)/sizex));
f2=0.5*(1+cos(2*pi*(t2-sizey/2-1)/sizey));
maskf=f1'*f2;
mask=zeros(dim(1),dim(2));
mask((dim(1)/2+1-sizex/2):(dim(1)/2+sizex/2),(dim(2)/2+1-sizey/2):(dim(2)/2+sizey/2))=maskf;
if length(dim)==2
    signal=mag.*exp(1i*phase);
    ks=fftshift(fftn(signal));
    k0=ks(dim(1)/2+1,dim(2)/2+1);
    ks=ks.*mask;
    ks(dim(1)/2+1,dim(2)/2+1)=k0;
    signaln=ifftn(fftshift(ks));
    hpima=angle(signal./signaln);
else
    h=waitbar(0,'Homodyne high-pass filtering');
    for i=1:dim(3)
        tempmag=mag(:,:,i);
        tempph=phase(:,:,i);
        signal=tempmag.*exp(1i*tempph);
        ks=fftshift(fftn(signal));
        ks=ks.*mask;
        signaln=ifftn(fftshift(ks));
        hpima(:,:,i)=angle(signal./signaln);
        SWIM_waitbar(i/dim(3),h);
    end
    close(h);
end
end