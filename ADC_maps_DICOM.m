% Copyright (C) Siamak P. Nejad-Davarani, University of Miami - All Rights Reserved
% You may use, distribute and modify this code under the terms of the CC BY-NC-SA 4.0 license.
% 
% 
% 
% The CC BY-NC-SA 4.0 license may be reviewed at https://creativecommons.org/licenses/by-nc-sa/4.0/


%Calculating ADC maps using all the selected directories

clear, clc
close all

numslice = 12;

%oldpath = cd;
%selpath = uigetdir;
allpaths = uigetfile_n_dir;

numdir = length(allpaths);

for f = 1:numdir

    selpath = char(allpaths{f})
    f
    selpath

if selpath==0
    error('Choose the folder!')
    return
end
%cd(selpath)




% List DICOM files in directory
dicom_files = rdir(strcat(selpath, '\**\*.IMA'));
% dicom_files = dir(strcat(selpath, '\*.dcm'));

% List DICOM files in directory
%dicom_files = dir('*.dcm');
%dicom_files = {dicom_files.name};
b_buffer=[];

%creating the 4D image matrix
timage = dicomread(strcat(selpath, '\', dicom_files(1).name));
%timage = dicomread(dicom_files(1).name);
imsize = size(timage);
numfiles = length(dicom_files);
All_im = zeros(imsize(1), imsize(2), numslice, numfiles/numslice);


count = 1;
for i=1:numfiles/numslice
    for j = 1:numslice
    metadata = dicominfo(strcat(selpath, '\', dicom_files(count).name));
    b_buffer(metadata.InstanceNumber)=[metadata.Private_0019_100c];
    slnum = mod(metadata.InstanceNumber, numslice);
    if mod(metadata.InstanceNumber, numslice) == 0
        slnum = numslice;
    end

    All_im(:,:,slnum,ceil(metadata.InstanceNumber/numslice)) = dicomread(strcat(selpath, '\', dicom_files(count).name));
    count = count+1;
    end
end
b_buffer;
b_buffer_fin = b_buffer(1:numslice:end)
%b_buffer=unique(b_buffer)

%This was added to normalize eaxh group to itself
All_im(:,:,:,1:end) = All_im(:,:,:,1:end)./All_im(:,:,:,1);

%if (f==1)
%All_im_fin = All_im;
%b_buffer_fin = b_buffer;
%end

%if (f>1)
%All_im_fin = cat(4,All_im_fin,All_im);
%b_buffer_fin = cat(2,b_buffer_fin,b_buffer);
%end

All_im_fin = All_im;


end

Num_b=length(b_buffer_fin);
ADC = zeros(imsize(1), imsize(2), numslice);

for x = 1:imsize(1)
    for y = 1:imsize(2)
        for z = 1:numslice
   
        sig = squeeze(All_im_fin(x,y,z,:));
        %plot(b_buffer,sig)
        %p = polyfit(-b_buffer_fin,log(sig/sig(1)),1);
        %figure(1), plot(b_buffer_fin, sig/sig(1), '*')
        p = polyfit(-b_buffer_fin,log(sig),1);
        %figure(1), plot(b_buffer_fin, sig, '*')

        ADC(x,y,z) = 1000*p(1);
        %pause
        
        end
    end
end

fname = strcat(selpath, '');
fid = fopen(fname, 'w');
fwrite(fid, ADC, 'float32')
fclose(fid)