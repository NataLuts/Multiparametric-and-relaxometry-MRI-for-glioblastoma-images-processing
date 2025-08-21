function [rs,ri]=GetScaleParadi(di)
rs=1;
ri=0;
if isfield(di,'RescaleSlope')
    rs=di.RescaleSlope;
end
if isfield(di,'RescaleIntercept')
    ri=di.RescaleIntercept;
end 