% Concatenate the files

File1 = 'E:\Images\ViewRay\STAGE\New11212018\VR_BRAIN_003_WEEK0_SIM\Analyze\04-Unknownx1.img'
File2 = 'E:\Images\ViewRay\STAGE\New11212018\VR_BRAIN_003_WEEK0_SIM\Analyze\04-Unknownx2.img'
File3 = 'E:\Images\ViewRay\STAGE\New11212018\VR_BRAIN_003_WEEK0_SIM\Analyze\04-Unknownx3.img'
File4 = 'E:\Images\ViewRay\STAGE\New11212018\VR_BRAIN_003_WEEK0_SIM\Analyze\All_T1.img'


fin = zeros(208, 256, 64*3);

fid = fopen(File1, 'r','l');
count =1;
for i = 1:64
    
A = fread(fid, [208, 256], 'uint16');
fin(:,:,count) = A;
max(max(A))
imagesc(A), colormap(gray)

pause
count = count+1;    
end

fclose(fid)


fid = fopen(File2, 'r','l');
for i = 1:64
    
A = fread(fid, [208, 256], 'uint16');
fin(:,:,count) = A;
max(max(A))
imagesc(A), colormap(gray)

pause
count = count+1;    
end

fclose(fid)

fid = fopen(File3, 'r','l');
for i = 1:64
    
A = fread(fid, [208, 256], 'uint16');
fin(:,:,count) = A;
max(max(A))
imagesc(A), colormap(gray)

pause
count = count+1;    
end


fid = fopen(File4, 'w','l');

fwrite(fid, fin, 'uint16');

fclose(fid)


