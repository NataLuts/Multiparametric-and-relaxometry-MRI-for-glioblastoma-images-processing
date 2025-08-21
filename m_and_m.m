function result = m_and_m(dataset, overslices, slidingwindow, start, numread, option)
%if option==1, mIP; else MIP.
    dim = size(dataset);
    if slidingwindow==0
       slidingwindow=overslices;
    end   
    ind_start=start:slidingwindow:(start+numread-1);
    result=zeros(dim(1),dim(2),length(ind_start),'single');
    for i=1:length(ind_start)
       if ind_start(i)+overslices-1>dim(3)
           ind_end=dim(3);
       else
           ind_end=ind_start(i)+overslices-1;
       end
       if option==1
           result(:,:,i)=min(dataset(:,:,ind_start(i):ind_end),[],3);
       else
           result(:,:,i)=max(dataset(:,:,ind_start(i):ind_end),[],3);
       end
    end
end
