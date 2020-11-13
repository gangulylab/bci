function [fp,fn,tp,tn] = balan_acc(p,grp_test)
% function [fp,fn,tp,tn] = balan_acc(p,grp_test)

% balanced accuracy
tp=0;tn=0;fp=0;fn=0;
for ii=1:length(p)
    if p(ii)==1 && grp_test(ii)==1
        tp=tp+1;
    end
    
    if p(ii)==-1 && grp_test(ii)==-1
        tn=tn+1;
    end
    
    if p(ii)==-1 && grp_test(ii)==1
        fn=fn+1;
    end
    
    if p(ii)==1 && grp_test(ii)==-1
        fp=fp+1;
    end
    
end
