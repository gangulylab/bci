function out = pool_features(temp,chmap)

% get the pooled data
new_temp=[];
[xx yy] = size(chmap);
for k=1:size(temp,2)
    tmp1 = temp(1:128,k);tmp1 = tmp1(chmap);
    tmp2 = temp(129:256,k);tmp2 = tmp2(chmap);
    tmp3 = temp(257:384,k);tmp3 = tmp3(chmap);
    pooled_data=[];
    for i=1:2:xx
        for j=1:2:yy
            delta = (tmp1(i:i+1,j:j+1));delta=mean(delta(:));
            beta = (tmp2(i:i+1,j:j+1));beta=mean(beta(:));
            hg = (tmp3(i:i+1,j:j+1));hg=mean(hg(:));
            pooled_data = [pooled_data; delta; beta ;hg];
        end
    end
    new_temp= [new_temp pooled_data];
end
out=new_temp;