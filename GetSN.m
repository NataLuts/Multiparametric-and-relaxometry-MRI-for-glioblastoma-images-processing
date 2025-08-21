function strSN = GetSN(filename)

idx=strfind(filename,'-');
maxidx=max(idx);
if maxidx==0
    strSN=filename;
else
    strSN=filename(1:maxidx);
end