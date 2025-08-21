function newimg = read_analyze(filepathname)
% useage: read_analyze('full file path including the extension (of the .hdr file)')
try
    nii=load_nii(filepathname);
catch ME
    nii=load_untouch_nii(filepathname);
end
img=nii.img;
dim=size(img);
if length(dim)==3
    newimg=zeros(dim(2),dim(1),dim(3));
    for i=1:dim(3)
        newimg(:,:,i)=flipud(img(:,:,i)');
    end
elseif length(dim)==2
    newimg=flipud(img');
elseif length(dim)==4
    newimg=zeros(dim(2),dim(1),dim(3),dim(4));
    for ii=1:dim(4)
        for i=1:dim(3)
            newimg(:,:,i,ii)=flipud(img(:,:,i,ii)');
        end
    end
end