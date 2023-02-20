%% SMOOTH BATCH KALMAN FILTER CODE
% Step 1: Train an KF's parameters from imagined movement data
% Step 2: Update the KF's parameters from closed-loop data with ReFit
%         assumption

%%%%%%% KALMAN FILTER %%%
% equations for KF
% model: Xhat(t+1) = A*Xhat(t) + K[y(t)-C*A*Xhat(t)]
% state: X(t+1) = A*X(t) + w(t) ; need a W noise cov matrix
% measurements: y(t) = C*X(t) + q(t); need a Q noise cov matrix
% C, K and Q estimated 
% W is set a priori
% A is also set a priori

%% STEP 1: TRAINING KF FROM IMAGINED MOVEMENT DATA

clc;
clear;
close all

%%%% CHANGE THIS ACCORDINGLY 
folderpath = 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20220520\Robot';
folders={'141735','142458','143022','144223','144607','145327','145327'};
%%%%%


load(fullfile('clicker','kf_robot_params.mat'));


files=[];
for i=1:length(folders)
    filepath = fullfile(folderpath,folders{i},'Imagined');
    files=[files;findfiles('',filepath)'];
end

% get the training data
neural=[];
kinematics=[];
for i=1:length(files)

    % load the data
    disp(i/length(files)*100)
    load(files{i})
    idx=find(TrialData.TaskState==3);
    kin = TrialData.CursorState;
    kin = kin(:,idx);
    neural_features = TrialData.SmoothedNeuralFeatures;
    neural_features = cell2mat(neural_features(idx));
    
    % get only non-zero velocities
    idx=abs(sum(kin(4:6,:),1))>0;
    kin = kin(:,idx);
    neural_features = neural_features(:,idx);    
    
    % delta, beta and hG
    fidx = [129:256 513:640 769:896]; 
    neural_features = neural_features(fidx,:);    

    % pooling
    neural_features = pool_features(neural_features,TrialData.Params.ChMap);  

    % store
    neural = [neural neural_features];
    kinematics = [kinematics kin];     
end

% STATES
X=kinematics;
X(7,:)=1;

% MEASUREMENTS
Y=neural;

% A matrix
A = KF_robot.A;

% W matrix
W = KF_robot.W;

% least squares estimates of velocity C
C = Y/X(4:end,:); 
C = [zeros(size(C,1),3) C];

% get estimates of Q 
q = Y-C*X;
Q = cov(q');
%Q = (1/size(Y,2)-1)* (q*q');

% compute Kalman Gain and Error Covariance 
P = rand(7);
K = randn(size(C'));
chk=1;chk1=1;
counter=0;
norm_val=[];
norm_val1=[];
iter=1;
while (chk>1e-10 && chk1>1e-10) || (iter<1000)
    temp_norm = norm(P(:));  
    temp_norm1 = norm(K(:));  
    P = A*P*A' + W;
    P(1:3,:) = zeros(3,size(P,2));
    P(:,1:3) = zeros(size(P,1),3);
    P(:,end)=0;
    P(end,:)=0;
    K = P*C'*pinv(C*P*C' + Q);
    P = P - K*C*P;
    chk = abs(temp_norm - norm(P(:)));
    chk1 = abs(temp_norm1 - norm(K(:)));
    counter=counter+1;
    norm_val = [norm_val chk];
    norm_val1 = [norm_val1 chk1];
    iter=iter+1;
end

% store results
KF_robot.C = C;
KF_robot.Q = Q;
KF_robot.P = P;
KF_robot.K = K;

save(fullfile('clicker','kf_robot_params'),'KF_robot','-v7.3','-nocompression');


%% STEP 2: UPDATING KF FROM CLOSED-LOOP ONLINE DATA USING SMOOTH BATCH

clc;
clear;
close all

folderpath = 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20220520\Robot';
folders={'141735','142458','143022','144223','144607','145327','145327'};

load(fullfile('clicker','kf_robot_params.mat'));


files=[];
for i=1:length(folders)
    filepath = fullfile(folderpath,folders{i},'BCI_Fixed');
    files=[files;findfiles('',filepath)'];
end

% ReFit: rotate decoded velocity vector towards target (intended vector)
% get the training data
neural=[];
kinematics=[];
for i=1:length(files)
    
    % load the data
    disp(i/length(files)*100)
    load(files{i})
    idx=find(TrialData.TaskState==3);
    kin = TrialData.CursorState;
    kin = kin(:,idx);
    neural_features = TrialData.SmoothedNeuralFeatures;
    neural_features = cell2mat(neural_features(idx));
    
    % get only non-zero velocities
    idx=abs(sum(kin(4:6,:),1))>0;
    kin = kin(:,idx);
    neural_features = neural_features(:,idx);    

    % rotate decoded velocity vector towards target
    target_pos = TrialData.TargetPosition;
    xint = [];
    K = TrialData.Params.KF_robot.K;
    C = TrialData.Params.KF_robot.C;
    for j=1:size(kin,2)

        % get current position
        pos = kin(1:3,j);
        recorded_vel = kin(4:6,j);
        % rotate towards target but keep same length as 
        new_vel_vector = target_pos' - pos;
        new_vel_vector = (new_vel_vector-)
       
       



    end
    
    % delta, beta and hG
    fidx = [129:256 513:640 769:896]; 
    neural_features = neural_features(fidx,:);    

    % pooling
    neural_features = pool_features(neural_features,TrialData.Params.ChMap);  

    % store
    neural = [neural neural_features];
    kinematics = [kinematics kin];     
end


% STATES
X=kinematics;
X(7,:)=1;

% MEASUREMENTS
Y=neural;

% A matrix
A = KF_robot.A;

% W matrix
W = KF_robot.W;

% least squares estimates of velocity C
C = Y/X(4:end,:); 
C = [zeros(size(C,1),3) C];

% get estimates of Q 
q = Y-C*X;
Q= cov(q');
%Q = (1/size(Y,2)-1)* (q*q');

% SMOOTH BATCH UPDATE OF C AND Q
alp = 0.75;
bet = 0.75;


% compute Kalman Gain and Error Covariance 
P = rand(7);
K = randn(size(C'));
chk=1;chk1=1;
counter=0;
norm_val=[];
norm_val1=[];
while chk>1e-10 && chk1>1e-10
    temp_norm = norm(P(:));  
    temp_norm1 = norm(K(:));  
    P = A*P*A' + W;
    P(1:3,:) = zeros(3,size(P,2));
    P(:,1:3) = zeros(size(P,1),3);
    P(:,end)=0;
    P(end,:)=0;
    K = P*C'*pinv(C*P*C' + Q);
    P = P - K*C*P;
    chk = abs(temp_norm - norm(P(:)));
    chk1 = abs(temp_norm1 - norm(K(:)));
    counter=counter+1;
    norm_val = [norm_val chk];
    norm_val1 = [norm_val1 chk1];
end

% store results
KF_robot.C = C;
KF_robot.Q = Q;
KF_robot.P = P;
KF_robot.K = K;

save(fullfile('clicker','kf_robot_params'),'KF_robot','-v7.3','-nocompression');




