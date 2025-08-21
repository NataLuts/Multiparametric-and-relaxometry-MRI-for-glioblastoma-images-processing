%Yongsheng Chen, andrew.ys.chen@gmail.com
%August 15, 2016
function [T1,p0,r]=fn_vfa_t1mapping(mag,FAs,k,TR,maxT1,maxPD,maxPDscale)
%least squre
%y=bx+a

dim=size(mag);
N=dim(4);
theta=FAs.*pi./180.0;
y=zeros(dim,'single');
x=zeros(dim,'single');

for i=1:dim(4)
   y(:,:,:,i)=mag(:,:,:,i)./sin(k(:,:,:,i).*theta(i));
   x(:,:,:,i)=mag(:,:,:,i)./tan(k(:,:,:,i).*theta(i));
end

sigma_x   = zeros(dim(1),dim(2),dim(3));
sigma_y   = zeros(dim(1),dim(2),dim(3));
sigma_xy  = zeros(dim(1),dim(2),dim(3));
sigma_sqx = zeros(dim(1),dim(2),dim(3));

for i=1:dim(4)
    sigma_x(:,:,:)   =  sigma_x(:,:,:)  +  x(:,:,:,i);
    sigma_y(:,:,:)   =  sigma_y(:,:,:)  +  y(:,:,:,i);
    sigma_xy(:,:,:)  =  sigma_xy(:,:,:) +  (x(:,:,:,i).*y(:,:,:,i));
    sigma_sqx(:,:,:) =  sigma_sqx(:,:,:)+  x(:,:,:,i).^2;
end

intercept = (sigma_sqx.*sigma_y-sigma_x.*sigma_xy)./((N.*sigma_sqx)-sigma_x.^2);
slope = (N.*sigma_xy-sigma_x.*sigma_y)./((N.*sigma_sqx)-sigma_x.^2);


mean_x = sigma_x/N;
mean_y = sigma_y/N;

sigma_avg=zeros(dim(1),dim(2),dim(3));
sigma_xpsq=zeros(dim(1),dim(2),dim(3));
sigma_ypsq=zeros(dim(1),dim(2),dim(3));

for i=1:dim(4)
    sigma_avg(:,:,:)=sigma_avg(:,:,:)+(x(:,:,:,i)-mean_x(:,:,:)).*(y(:,:,:,i)-mean_y(:,:,:));
    sigma_xpsq(:,:,:)=sigma_xpsq(:,:,:)+(x(:,:,:,i)-mean_x(:,:,:)).^2;
    sigma_ypsq(:,:,:)=sigma_ypsq(:,:,:)+(y(:,:,:,i)-mean_y(:,:,:)).^2;
end

r=sigma_avg./((sigma_xpsq.*sigma_ypsq).^0.5);


T1=-TR./log(slope);
T1=abs(T1);

p0=intercept./(1-slope);
p0=abs(p0);

p0=p0/maxPDscale;

T1(T1>maxT1)=maxT1;
T1(T1<0)=0;
T1(isnan(T1))=0;

p0(p0>maxPD)=maxPD;
p0(p0<0)=0;
p0(isnan(p0))=0;






