% Copyright (C) Siamak P. Nejad-Davarani, University of Miami - All Rights Reserved
% You may use, distribute and modify this code under the terms of the CC BY-NC-SA 4.0 license.
% 
% 
% 
% The CC BY-NC-SA 4.0 license may be reviewed at https://creativecommons.org/licenses/by-nc-sa/4.0/


%Calculating ADC maps using all the selected directories

clear, clc
close all

selpath = uigetdir

% This is the threshold that is used for eliminating the background 
z_thresh = 1;


%selpath = ''
%selpath = ''

% List DICOM files in directory
nifti_files = dir(strcat(selpath, '\*.hdr'))


timage = analyze75read(strcat(selpath, '\', nifti_files(1).name));
info = niftiinfo(strcat(selpath, '\', nifti_files(1).name));

imsize = size(timage);
numfiles = length(nifti_files);


All_im = zeros(imsize(1), imsize(2), imsize(3), numfiles);

bval_buffer = zeros(1,numfiles);


for i=1:numfiles
    fname = strcat(selpath, '\',nifti_files(i).name);
    b_buffer(i) = extractBV(fname);

    All_im(:,:,:, i) = analyze75read(strcat(selpath, '\',nifti_files(i).name));
   
   
end

b_buffer_fin = b_buffer
zers = find(b_buffer == 0)
zers(end+1) = length(b_buffer_fin)+1
All_im_f = zeros(size(All_im));

%normalization of the S matrices to each of the S0s in the series

for i = 1:length(zers)-1

    for j = zers(i):zers(i+1)-1 
     %[i, j, zers(i), zers(i+1)-1 ]
        s0 = All_im(:,:,:, zers(i));
        All_im_f(:,:,:, j) = All_im(:,:,:, j) ./s0;
    end

end



%pause
Num_b=length(b_buffer_fin);
ADC = zeros(imsize(1), imsize(2), imsize(3));

av_mat = mean(All_im,4);



for x = 1:imsize(1)
    x
    for y = 1:imsize(2)
        for z = 1:imsize(3)
   
            if av_mat(x,y,z) > z_thresh

        sig = squeeze(All_im_f(x,y,z,:));
        %figure(1), plot(b_buffer,sig/sig(1), '*')
        %p = polyfit(-b_buffer_fin,log(sig/sig(1)),1);
        %p = polyfit(-b_buffer_fin,log(sig),1);
        
        %Function #1 
                B = log(sig);
                A = -b_buffer_fin';

                AA = (A'*A);
                [U,S,V] = svd(AA,'econ');
                X1 = (V*inv(S)*U')*A'*B;
            
        
        %figure(1), plot(b_buffer,sig, '*')

        ADC(x,y,z) =X1;
        %pause
            end
        end
    end
end

fname = '';
fid = fopen(fname, 'w');
fwrite(fid, ADC, 'float32')
fclose(fid)
%% 
% vox_size = [1.5038 1.5038 1.5]
vox_size = [2.3438 2.3438 9]
Analyze = UT_Write_Analyze ('', ADC, 'float', vox_size);