function [rs,ri]=GetScalePara(strdcm)
rs=1;
ri=0;
di=dicominfo(strdcm);
if isfield(di,'RescaleSlope')
    rs=di.RescaleSlope;
end
if isfield(di,'RescaleIntercept')
    ri=di.RescaleIntercept;
end 