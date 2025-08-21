%Yongsheng Chen <andrew.ys.chen@gmail.com>
%July 16, 2018

clear;
clc;

%you may firstly use SPIN software to sort the dicom forlder exported from the
%scanner, or read them in your own way.
%dirPatient='C:\ProjectData\STAGE_CGH_VR_BRAIN_001_WEEK4\';
dirPatient='';


%dirA1='06-A1-PDW-mag\';
%dirA1p='07-A1-PDW-pha\';
dirA1='';
%dirA1p='09-Unknown\';

%dirA2='02-A2-T1W-mag\';
%dirA2p='03-A2-T1W-pha\';
dirA2='';
%dirA2p='05-Unknown\';


%dirTrufi='10-trueFISP-mag\';
dirTrufi='10-Unknown\';

%fileA1 = 'E:\Images\ViewRay\STAGE\VR_BRAIN_009_ALL\VR_BRAIN_009_WK2_coreg\rPD_all.img'
%fileA2 = 'E:\Images\ViewRay\STAGE\VR_BRAIN_009_ALL\VR_BRAIN_009_WK2_coreg\T1_all.img'
%fileTr = 'E:\Images\ViewRay\STAGE\VR_BRAIN_009_ALL\VR_BRAIN_009_WK2_coreg\r10-Unknown.img'

%usually I measure the noise level at the edge of the head to determin
%sNoise. Any region signal intensity less than sNoise will be discarded.
sNoise=35;

bSaveDicom = 1; % 1: to save dicom format; otherwise not to save dicom


% nMIP=4;

%% processing

dirA1=[dirPatient,dirA1];
%dirA1p=[dirPatient,dirA1p];

dirA2=[dirPatient,dirA2];
%dirA2p=[dirPatient,dirA2p];

%dirTrufi=[dirPatient,dirTrufi];

dirOutput=[dirPatient,'STAGE\'];
if 0==exist(dirOutput,'dir')
   mkdir(dirOutput); 
end

%dicom info
filesA1 = dir([dirA1, '\*.dcm']);
filesA2 = dir([dirA2, '\*.dcm']);
%filesTruefi = dir([dirTrufi, '\*.dcm']);

snA1 = GetSN(filesA1(1).name);
% snA1p = dirA1p(1:2);

snA2 = GetSN(filesA2(1).name);
% snA2p = dirA2p(1:2);

%snTrufi= GetSN(filesTruefi(1).name);

slice   = length(filesA1(not([filesA1.isdir])))/3; %three echoes

di1=dicominfo([dirA1,snA1,num2str(1+slice*0),'.dcm']);
di2=dicominfo([dirA1,snA1,num2str(1+slice*1),'.dcm']);
di3=dicominfo([dirA1,snA1,num2str(1+slice*2),'.dcm']);

di4=dicominfo([dirA2,snA2,num2str(1+slice*0),'.dcm']);

%di5=dicominfo([dirTrufi,snTrufi,'1.dcm']);

row=di1.Rows;
column=di1.Columns;
%row = 192;
%column = 256;
%slize = 64;

TE1=di1.EchoTime;
TE2=di2.EchoTime;
TE3=di3.EchoTime;

dim=[row,column, slice,3];
%dim = [256, 208, 64,3]
TEs=[TE1,TE2,TE3];

ddx=di1.PixelSpacing(1);
ddy=di1.PixelSpacing(2);
ddz=di1.SliceThickness;
ori=di1.ImageOrientationPatient;
iFreq=di1.ImagingFrequency;

faA1=di1.FlipAngle;
faA2=di4.FlipAngle;
TR=di1.RepetitionTime;

%swim parameters
% betThreshold=0.2;
% betErode=0;
% betIsland=2000;
% qmapThreshold=0.2;
% InverseFilterThreshold=0.1; 
% Iterations=4;
% IterationThreshold=0.15;
% KernalSize=6;
% uwType=2;

%orignial data
magA1=zeros(dim,'single');
phaseA1=zeros(dim,'single');

magA2=zeros(dim,'single');
phaseA2=zeros(dim,'single');

%magTrufi = zeros(row, column, slice);

%load data
[rs1,ri1]=GetScaleParadi(di1);
[rs2,ri2]=GetScaleParadi(di4);

%[rs3,ri3]=GetScaleParadi(di5);



%fidA1 = fopen(fileA1, 'r', 'b');
%fidA2 = fopen(fileA2, 'r', 'b');
%fidTr = fopen(fileTr, 'r', 'b');


for e=1:dim(4)
    for s=1:dim(3)
        magA1(:,:,s,e)=dicomread([dirA1,snA1,num2str(s+(dim(3)*(e-1))),'.dcm'])*rs1+ri1;
%         phaseA1(:,:,s,e)=dicomread([dirA1p,snA1p,'-',num2str(s+(dim(3)*(e-1))),'.dcm']);
        %magA1(:,:,s,e)= (fliplr(fread(fidA1, [dim(2), dim(1)], 'uint16'))')*rs1+ri1;
        %figure(1), imagesc(magA1(:,:,s,e)), colormap(gray)
        
        
        magA2(:,:,s,e)=dicomread([dirA2,snA2,num2str(s+(dim(3)*(e-1))),'.dcm'])*rs2+ri2;
%       phaseA2(:,:,s,e)=dicomread([dirA2p,snA2p,'-',num2str(s+(dim(3)*(e-1))),'.dcm']);
        %magA2(:,:,s,e)= (fliplr(fread(fidA2, [dim(2), dim(1)], 'uint16'))')*rs2+ri2;
        %figure(2), imagesc(magA2(:,:,s,e)), colormap(gray)
        
        if e==1
            %magTrufi(:,:,s)=dicomread([dirTrufi,snTrufi,num2str(s),'.dcm'])*rs3+ri3;
        %    magTrufi(:,:,s)=(fliplr(fread(fidTr, [dim(2), dim(1)], 'uint16'))')*rs3+ri3;
        %figure(3), imagesc(magTrufi(:,:,s,e)), colormap(gray)
        %pause
        
        end
    end
end

%fclose(fidA1)
%fclose(fidA2)
%fclose(fidTr)


%% TE2 swim, tswi, swi
% for now SWIM/QSM, SWI and tSWI are useless for 0.35T. therefor, I didn't
% include these processing for now.

% [swim_te2,iswim_te2,swi_te2,tswi_te2,betMask_te2,tswiMaks_te2,php_te2,maskn_te2] = fn_swim_swi_tswi(...
%     magA1(:,:,:,2),phaseA1(:,:,:,2),dim(1),dim(2),dim(3),TEs(2),ddx,ddy,ddz,ori,iFreq,...
%     betThreshold,betErode,betIsland,qmapThreshold,InverseFilterThreshold,Iterations,IterationThreshold,KernalSize,uwType);

% php_te2(magA1(:,:,:,1)<sNoise)=0;
% swim_mip=m_and_m(iswim_te2,nMIP,1,1,dim(3),2);
% pswim_mip=m_and_m(php_te2,nMIP,1,1,dim(3),2);
% swi_mip=m_and_m(swi_te2,nMIP,1,1,dim(3),1);
% tswi_mip=m_and_m(tswi_te2,nMIP,1,1,dim(3),1);

% mat2analyze_file(swi_te2,[dirOutput,'swi_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(tswi_te2,[dirOutput,'tswi_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(iswim_te2*1000,[dirOutput,'swim_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(php_te2/pi*4096,[dirOutput,'pswim_te2.hdr'],[ddx,ddy,ddz],16);
% 
% mat2analyze_file(swi_mip,[dirOutput,'mip4_swi_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(tswi_mip,[dirOutput,'mip4_tswi_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(swim_mip*1000,[dirOutput,'mip4_swim_te2.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(pswim_mip/pi*4096,[dirOutput,'mip4_pswim_te2.hdr'],[ddx,ddy,ddz],16);

%% TE3 swim, tswi, swi
% [swim_te3,iswim_te3,swi_te3,tswi_te3,betMask_te3,tswiMaks_te3,php_te3,maskn_te3] = fn_swim_swi_tswi(...
%     magA1(:,:,:,3),phaseA1(:,:,:,3),dim(1),dim(2),dim(3),TEs(3),ddx,ddy,ddz,ori,iFreq,...
%     betThreshold,betErode,betIsland,qmapThreshold,InverseFilterThreshold,Iterations,IterationThreshold,KernalSize,uwType);
% 
% php_te3(magA1(:,:,:,1)<sNoise)=0;
% swim_mip=m_and_m(iswim_te3,nMIP,1,1,dim(3),2);
% pswim_mip=m_and_m(php_te3,nMIP,1,1,dim(3),2);
% swi_mip=m_and_m(swi_te3,nMIP,1,1,dim(3),1);
% tswi_mip=m_and_m(tswi_te3,nMIP,1,1,dim(3),1);
% 
% mat2analyze_file(swi_te3,[dirOutput,'swi_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(tswi_te3,[dirOutput,'tswi_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(iswim_te3*1000,[dirOutput,'swim_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(php_te3/pi*4096,[dirOutput,'pswim_te3.hdr'],[ddx,ddy,ddz],16);
% 
% mat2analyze_file(swi_mip,[dirOutput,'mip4_swi_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(tswi_mip,[dirOutput,'mip4_tswi_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(swim_mip*1000,[dirOutput,'mip4_swim_te3.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(pswim_mip/pi*4096,[dirOutput,'mip4_pswim_te3.hdr'],[ddx,ddy,ddz],16);

%% t1 mapping
%t1 mapping, there is no RF penetration correction needed for 0.35T. 
%for high field B1+ and B1- correction, please refer to STAGE part-II paper

k=ones(row,column,slice,2);
mag=zeros(row,column,slice,2);
FAs=[faA1,faA2];
maxT1=20000;
maxPD=20000;
maxPDscale=1;

%TE1
mag(:,:,:,1)=magA1(:,:,:,1);
mag(:,:,:,2)=magA2(:,:,:,1);
[T1_te1,p0_te1,rsq_te1]=fn_vfa_t1mapping(mag,FAs,k,TR,maxT1,maxPD,maxPDscale);

%TE2
mag(:,:,:,1)=magA1(:,:,:,2);
mag(:,:,:,2)=magA2(:,:,:,2);
[T1_te2,p0_te2,rsq_te2]=fn_vfa_t1mapping(mag,FAs,k,TR,maxT1,maxPD,maxPDscale);

%average
T1_avg=(T1_te1+T1_te2)/2;
p0_avg=(p0_te1+p0_te2)/2;

T1_avg(magA1(:,:,:,1)<sNoise)=0;
p0_avg(magA1(:,:,:,1)<sNoise)=0;

mat2analyze_file(T1_avg,[dirOutput,'T1map.hdr'],[ddx,ddy,ddz],16);
mat2analyze_file(p0_avg,[dirOutput,'PDmap.hdr'],[ddx,ddy,ddz],16);

%% pure PDMAP
% deltaE2=p0_te1./p0_te2;
% deltaE2(isnan(deltaE2))=0;
% deltaE2(deltaE2>1000)=0;
% 
% E2TE1=deltaE2.^(-1*TE1/(TE2-TE1));
% 
% purePD=p0_te1./E2TE1;
% purePD(purePD<1)=0;
% purePD(isnan(purePD))=0;
% purePD(purePD>20000)=20000;
% mat2analyze_file(purePD,[dirOutput,'purePD.hdr'],[ddx,ddy,ddz],16);

%% R2star map
%the delta TE (TE3-TE1) is too short for 0.35T, therefore, the r2star and
%t2star maps should have very low SNR.

[R2sA1, protonDensityA1,maskA1] = STAGE_R2star(magA1, TEs, 3, ones(size(magA1(:,:,:,1))));
[R2sA2, protonDensityA2,maskA2] = STAGE_R2star(magA2, TEs, 3, ones(size(magA2(:,:,:,1))));

R2sA1(magA1(:,:,:,1)<sNoise)=-1;
R2sA2(magA1(:,:,:,1)<sNoise)=-1;

avgR2star=(R2sA2+R2sA1)/2;
% mat2analyze_file(R2sA1*1000,[dirOutput,'r2s_A1.hdr'],[ddx,ddy,ddz],16);
% mat2analyze_file(R2sA2*1000,[dirOutput,'r2s_A2.hdr'],[ddx,ddy,ddz],16);
mat2analyze_file(avgR2star*1000,[dirOutput,'r2s_avg.hdr'],[ddx,ddy,ddz],16);

T2star=1/avgR2star;
T2star(T2star>200)=200;
T2star(T2star<1)=0;
mat2analyze_file(T2star,[dirOutput,'t2s_avg.hdr'],[ddx,ddy,ddz],16);

%% T1WE

T1WE1=magA2(:,:,:,1)-magA1(:,:,:,1);
T1WE2=magA2(:,:,:,2)-magA1(:,:,:,2);
T1WE3=magA2(:,:,:,3)-magA1(:,:,:,3);

T1WE=(T1WE1+T1WE2+T1WE3)/2;
minT1WE=abs(min(min(min(T1WE))));
T1WE=T1WE+minT1WE;
T1WE(magA1(:,:,:,1)<sNoise) = 0;

mat2analyze_file(T1WE,[dirOutput,'T1WE.hdr'],[ddx,ddy,ddz],16);

%% pT2W
csf=(magA1(:,:,:,1) + magA1(:,:,:,2) + magA1(:,:,:,3) + magA2(:,:,:,1) + magA2(:,:,:,2) + magA2(:,:,:,3))/6;

maxcsf=max(max(max(csf)));
mincsf=min(min(min(csf)));
pT2W=(1-(csf-mincsf)/(maxcsf-mincsf))*1000;
pT2W(magA1(:,:,:,1)<sNoise) = 500;

mat2analyze_file(pT2W,[dirOutput,'pT2W.hdr'],[ddx,ddy,ddz],16);


%% T2 mapping
%This T2 mapping would be Niloufar's current project. For now, I use the
%ratio of SSFP and T1 to get an assumed T2 scale.

%minp=min(min(min(magTrufi)));
%maxp=max(max(max(magTrufi)));
%norTrufi = (magTrufi - minp)/(maxp - minp);

% minp=min(min(min(p0_avg)));
% maxp=max(max(max(p0_avg)));
% norPD = (p0_avg - minp)/(maxp - minp);

% T2map=4*magTrufi.^2.*T1_avg./(norPD.^2);
%T2map=T1_avg.*norTrufi;
%mat2analyze_file(T2map,[dirOutput,'T2map.hdr'],[ddx,ddy,ddz],16);

%% simuFLAIR

TRir=10000;%ms
TEir=15;%ms 
TIir=2000;%ms

t1f=1-exp(-1);
t2f=exp(-1);

thNeg = -8191;
thPos = 8192;

T1forir    = T1_avg;
p0forir    = p0_avg;
%T2forir    = T2map;

%flair = p0forir .* (1-2*exp(-TIir./T1forir)) .* (1-exp(-TRir./T1forir)) .* (exp(-TEir./T2forir));
%flair(flair > thPos) = thPos;
%flair(flair < thNeg) = thNeg;
%mat2analyze_file(flair,[dirOutput,'simFLAIR.hdr'],[ddx,ddy,ddz],16);

%% Saving dicom files

if bSaveDicom == 1
    snSTAGE_T1WE=8001;
    snSTAGE_T1MAP=8002;
    snSTAGE_PDMAP=8003;
    snSTAGE_R2star=8004;
    snSTAGE_T2star=8005;
    snSTAGE_pT2W=8006;
    snSTAGE_simFLAIR=8007;
    snSTAGE_T2MAP=8008;

    dsSTAGE_T1WE='STAGE enhanced T1W';
    dsSTAGE_T1MAP='STAGE T1 MAP';
    dsSTAGE_PDMAP='STAGE PD MAP';
    dsSTAGE_R2star='STAGE R2 star MAP';
    dsSTAGE_T2star='STAGE T2 star MAP';
    dsSTAGE_pT2W='STAGE pseudo T2W';
    dsSTAGE_simFLAIR='STAGE simulated FLAIR';
    dsSTAGE_T2MAP='STAGE T2 MAP';

    dirSTAGE_T1WE=[dirOutput,'01_T1WE\'];
    dirSTAGE_T1MAP=[dirOutput,'02_T1MAP\'];
    dirSTAGE_PDMAP=[dirOutput,'03_PDMAP\'];
    dirSTAGE_R2star=[dirOutput,'04_R2star\'];
    dirSTAGE_T2star=[dirOutput,'05_T2star\'];
    dirSTAGE_pT2W=[dirOutput,'06_pT2W\'];
    dirSTAGE_simFLAIR=[dirOutput,'07_simFLAIR\'];
    dirSTAGE_T2MAP=[dirOutput,'08_T2MAP\'];

    if 0 == exist(dirSTAGE_T1WE,'dir') 
        mkdir(dirSTAGE_T1WE);
    end

    if 0 == exist(dirSTAGE_T1MAP,'dir')
        mkdir(dirSTAGE_T1MAP);
    end

    if 0 == exist(dirSTAGE_PDMAP,'dir')
        mkdir(dirSTAGE_PDMAP);
    end

    if 0 == exist(dirSTAGE_pT2W,'dir')
       mkdir(dirSTAGE_pT2W); 
    end

    if 0 == exist(dirSTAGE_simFLAIR,'dir')
     %  mkdir(dirSTAGE_simFLAIR); 
    end

    if 0 == exist(dirSTAGE_R2star,'dir') 
        mkdir(dirSTAGE_R2star);
    end

    if 0 == exist(dirSTAGE_T2star,'dir') 
        mkdir(dirSTAGE_T2star);
    end

    if 0 == exist(dirSTAGE_T2MAP,'dir')
     %   mkdir(dirSTAGE_T2MAP);
    end

    [ww_t1we,wc_t1we]=GetWindowLevel(T1WE);
    [ww_t1,wc_t1]=GetWindowLevel(T1_avg);
    [ww_p0,wc_p0]=GetWindowLevel(p0_avg);
    [ww_pt2w,wc_pt2w]=GetWindowLevel(pT2W);
    %[ww_simflair,wc_simflair]=GetWindowLevel(flair);
    [ww_r2s,wc_r2s]=GetWindowLevel(avgR2star*1000);
    [ww_t2s,wc_t2s]=GetWindowLevel(T2star);
    %[ww_t2,wc_t2]=GetWindowLevel(T2map);

    h1 = waitbar(0);
    strDateTime=datestr(now,'yyyymmddHHMMSSFFF');
    dicomRootUID='1.3.12.2.1107.5.2.36.40292.';% WSU MRRF Verio Scanner

    for i=1:slice
        waitbar((i) / slice, h1, 'writing dicom file......');

        di=dicominfo([dirA1,snA1,num2str(i),'.dcm']);
        
        %t1we
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_T1WE),'.',strDateTime];
        di.SeriesNumber=snSTAGE_T1WE;
        di.SeriesDescription=dsSTAGE_T1WE;
        di.ProtocolName=dsSTAGE_T1WE;
        di.WindowCenter=wc_t1we;
        di.WindowWidth=ww_t1we;
        dicomwrite(int16(T1WE(:,:,i)),strcat(dirSTAGE_T1WE,'T1WE-',num2str(i),'.dcm'),di);

        %t1 map
        di=dicominfo([dirA1,snA1,num2str(i),'.dcm']);
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_T1MAP),'.',strDateTime];
        di.SeriesNumber=snSTAGE_T1MAP;
        di.SeriesDescription=dsSTAGE_T1MAP;
        di.ProtocolName=dsSTAGE_T1MAP;
        di.WindowCenter=wc_t1;
        di.WindowWidth=ww_t1;
        dicomwrite(int16(T1_avg(:,:,i)),strcat(dirSTAGE_T1MAP,'T1MAP-',num2str(i),'.dcm'),di);

        %pd map
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_PDMAP),'.',strDateTime];
        di.SeriesNumber=snSTAGE_PDMAP;
        di.SeriesDescription=dsSTAGE_PDMAP;
        di.ProtocolName=dsSTAGE_PDMAP;
        di.WindowCenter=wc_p0;
        di.WindowWidth=ww_p0;
        dicomwrite(int16(p0_avg(:,:,i)),strcat(dirSTAGE_PDMAP,'PDMAP-',num2str(i),'.dcm'),di);
        
        %pseudo t2w
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_pT2W),'.',strDateTime];
        di.SeriesNumber=snSTAGE_pT2W;
        di.SeriesDescription=dsSTAGE_pT2W;
        di.ProtocolName=dsSTAGE_pT2W;
        di.WindowCenter=wc_pt2w;
        di.WindowWidth=ww_pt2w;
        dicomwrite(int16(pT2W(:,:,i)),strcat(dirSTAGE_pT2W,'pT2W-',num2str(i),'.dcm'),di);

        %simulated flair
   %     di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_simFLAIR),'.',strDateTime];
   %     di.SeriesNumber=snSTAGE_simFLAIR;
   %     di.SeriesDescription=dsSTAGE_simFLAIR;
   %     di.ProtocolName=dsSTAGE_simFLAIR;
   %     di.WindowCenter=wc_simflair;
   %     di.WindowWidth=ww_simflair;
   %     dicomwrite(int16(flair(:,:,i)),strcat(dirSTAGE_simFLAIR,'simFLAIR-',num2str(i),'.dcm'),di);

        %r2 star mpa
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_R2star),'.',strDateTime];
        di.SeriesNumber=snSTAGE_R2star;
        di.SeriesDescription=dsSTAGE_R2star;
        di.ProtocolName=dsSTAGE_R2star;
        di.WindowCenter=wc_r2s;
        di.WindowWidth=ww_r2s;
        dicomwrite(int16(avgR2star(:,:,i).*1000),[dirSTAGE_R2star,'R2star-',num2str(i),'.dcm'],di);

        %t2 star map
        di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_T2star),'.',strDateTime];
        di.SeriesNumber=snSTAGE_T2star;
        di.SeriesDescription=dsSTAGE_T2star;
        di.ProtocolName=dsSTAGE_T2star;
        di.WindowCenter=wc_t2s;
        di.WindowWidth=ww_t2s;
        dicomwrite(int16(T2star(:,:,i)),[dirSTAGE_T2star,'T2star-',num2str(i),'.dcm'],di);

        %t2 map
 %       di.SeriesInstanceUID=[dicomRootUID,num2str(snSTAGE_T2MAP),'.',strDateTime];
 %       di.SeriesNumber=snSTAGE_T2MAP;
 %       di.SeriesDescription=dsSTAGE_T2MAP;
 %       di.ProtocolName=dsSTAGE_T2MAP;
 %       di.WindowCenter=wc_t2;
 %       di.WindowWidth=ww_t2;
 %       dicomwrite(int16(T2map(:,:,i)),[dirSTAGE_T2MAP,'T2MAP-',num2str(i),'.dcm'],di);
    end

    close (h1);
end
