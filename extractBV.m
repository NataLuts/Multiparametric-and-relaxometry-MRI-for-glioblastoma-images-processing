% Copyright (C) Siamak P. Nejad-Davarani, University of Miami - All Rights Reserved
function te = extractBV(fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fname
    k = strfind(fname,"__");
    k2 = strfind(fname, "__bvalue");
    a1 = find(k==k2);
    ss = extractAfter(fname, k(a1)+1);
    k = strfind(ss,"__");
    ss = extractBefore(ss, k(1))
    f = strfind(ss, 'bvalue');
    te = str2num(extractAfter(ss,f+5));

end