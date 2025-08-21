function maskn=imdilate3d(mask,radius,dx,dy,dz)
%3D dilation of the mask
dim=size(mask);
mask2=zeros(dim(1),dim(2),dim(3)+2,'single');
mask2(:,:,2:dim(3)+1)=mask;
mask=mask2;
clear mask2;
dim=size(mask);
rho=spkernel(dim,radius,dx,dy,dz,1);
rhoft=fftn(rho);
maskn=fftshift(ifftn(fftn(mask).*rhoft));
maskn=round(maskn*10000)/10000;
maskn(maskn>0)=1;
maskn=maskn(:,:,2:dim(3)-1);
clear rho_sp rho rhoft
end