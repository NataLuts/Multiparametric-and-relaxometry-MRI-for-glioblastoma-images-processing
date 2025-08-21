function [w,c] = GetWindowLevel(data)
    minvalue=min(min(min(data)));
    maxvalue=max(max(max(data)));
    
    w = maxvalue-minvalue;
    c = (maxvalue+minvalue)/2;
    
end