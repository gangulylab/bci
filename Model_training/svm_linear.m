function [res_acc,model_overall,pval,model_runs] = svm_linear(A,B,iterations,prop)
%function [res_acc,model_overall,pval] = svm_linear(A,B,iterations,prop)
% rows - features
% columns - time-points or observations 


% is there a difference b/w the two cognitive states?
% via SVM
dd = pwd;
imagine_data=A;
non_imagine_data=B;
% linear SVM
cd('C:\Users\Nikhlesh\Documents\liblinear-2.1\liblinear-2.1\matlab')
% preparing stuff
% make the smaller dataset to be first i.e. imagine data and the larger
% dataset to be non_imagine_data
s1 = size(imagine_data,2);
s2 = size(non_imagine_data,2);


if s1<= s2
    m=s1;
else
    m=s2;
end

% getting trials together and group membership
if s1<=s2
    trials = [imagine_data';non_imagine_data'];
    grp = [-1*ones(size(imagine_data,2),1);ones(size(non_imagine_data,2),1)];
elseif s1>s2
    trials = [non_imagine_data';imagine_data'];
    grp = [-1*ones(size(non_imagine_data,2),1);ones(size(imagine_data,2),1)];
end

%%%%scaling SVM
%minn=min(trials);
%trials=trials-repmat(minn,[size(trials,1) 1]);
%maxx=max(trials);
%trials = trials./ repmat(maxx,[size(trials,1) 1]);

%%%% variables
model_overall=[];res_acc=[];
model_runs={};
parfor iter=1:iterations
    
    % randomly partitioning the data into training and test trials
    temp=round(prop*m);
    grpi_trial_index = randperm(m,temp);
    if s2<=s1
        grpni_trial_index = s2+randperm(s1,temp);
    else
        grpni_trial_index = s1+randperm(s2,temp);
    end
    grpi_trial=trials(grpi_trial_index,:);
    grpni_trial=trials(grpni_trial_index,:);
    train_trials=[grpi_trial;grpni_trial];
    grp_train=[-1*ones(size(grpi_trial,1),1); ones(size(grpni_trial,1),1)];
    
    %separating training trials from overall dataset
    temp1=[grpi_trial_index grpni_trial_index]';
    temp=ones(length(grp),1);
    temp(temp1)=0;
    temp=logical(temp);
    grp_test = grp(temp);
    test_trials=trials(temp,:);
    
    % scaling SVM
    %     minn=min(train_trials);
    %     train_trials=train_trials-repmat(minn,[size(train_trials,1) 1]);
    %     maxx=max(train_trials);
    %     train_trials = train_trials./ repmat(maxx,[size(train_trials,1) 1]);
    %     test_trials=test_trials-repmat(minn,[size(test_trials,1) 1]);
    %     test_trials = test_trials./ repmat(maxx,[size(test_trials,1) 1]);
    
    %training the SVM
    %to get the C parameter via simple grid search and train the model on
    %entire training set
    % C=[0.01:0.01:1 1:2:100];acc=[];
    % for j=1:length(C)
    %     par=['-c ' num2str(C(j)) ' -s 2 -v 10 -q']; % SVM parameters
    %     acc=[ acc train(grp_train,sparse(train_trials),par)];
    % end
    
    % find C using inbuilt function
    par=['-C -s 2 -v 10 -q']; % SVM parameters
    res=train(grp_train,sparse(train_trials),par);
    
    %[t1 t2]=max(acc);
    %par=['-c ' num2str(C(t2)) ' -s 2 '];
    par=['-c ' num2str(res(1)) ' -s 2 '];
    model=train(grp_train,sparse(train_trials),par);
    
    % test the model
    [p, a ,d]=predict(grp_test,sparse(test_trials),model);
    model_runs(iter).p=p;
    model_runs(iter).d=d;
    model_runs(iter).grp_test=grp_test;
    model_runs(iter).test_trials=test_trials;
    
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
    res_acc(iter) = 0.5* ( tp/(tp+fn) + tn/(tn+fp) );
    if s1<=s2
        model_overall(iter,:) = model.w;    
    else
        model_overall(iter,:) = -model.w;    
    end
    
    % p-vaue by bayesian stats    
    alp1=1+tp;
    bet1=1+fn;
    alp2=1+tn;
    bet2=1+fp;
    res=0.001;
    u=0:res:1;
    a=betapdf(u,alp1,bet1);
    b=betapdf(u,alp2,bet2);
    x=conv(a,b);
    z=2*x(1:2:end);
    z=z/(sum(x*res));
   % figure;plot(u,z);hold on;plot(u,a,'k');plot(u,b,'r')
    % calculate p-value
    querypt= 0.5;
    I=(u>querypt);
    pval(iter)=1-sum(z(I)*res);
    
     model_runs(iter).z=z;
     model_runs(iter).u=u;
     model_runs(iter).a=a;
     model_runs(iter).b=b;
    
end
cd(dd)

%mean(res_acc)*100
%temp=model_overall;
%figure;stem(mean(temp,1))