%% SMOOTH BATCH KALMAN FILTER CODE
% Step 1: Train an KF's parameters from imagined movement data
% Step 2: Update the KF's parameters from closed-loop data with ReFit
%         assumption

%%%%%%% KALMAN FILTER %%%
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

addpath(fullfile('task_helpers'))

load(fullfile('clicker','kf_robot_params.mat'));

ch_pooling=false;

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
    if ch_pooling == true
        neural_features = pool_features(neural_features,TrialData.Params.ChMap);  
    end

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
while chk>1e-10 || chk1>1e-10 
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
    if iter==1000
        break
    end
end

% store results
KF_robot.C = C;
KF_robot.Q = Q;
KF_robot.P = P;
KF_robot.K = K;

save(fullfile('clicker','kf_robot_params'),'KF_robot','-v7.3','-nocompression');
disp('Kalman Filter seeded on imagined data')
clearvars -except KF_robot

%% STEP 2: UPDATING KF FROM CLOSED-LOOP ONLINE DATA USING SMOOTH BATCH

clc;
clear;
close all

addpath(fullfile('task_helpers'))

%%%%% CHANGE AS REQUIRED %%%%%
folderpath = 'F:\DATA\ecog data\ECoG BCI\GangulyServer\Multistate clicker\20220520\Robot';
folders={'141735','142458','143022','144223','144607','145327','145327'};
%%%%% CHANGE AS REQUIRED %%%%%


load(fullfile('clicker','kf_robot_params.mat'));


files=[];
for i=1:length(folders)
    filepath = fullfile(folderpath,folders{i},'BCI_Fixed');
    files=[files;findfiles('',filepath)'];
end


% pooling option
ch_pooling=false;

% SMOOTH BATCH UPDATES: BETWEEN 0 AND 1 -> HOW MUCH THE PREVIOUS KF SHOULD
% BE WEIGHTED AS COMPARED TO ONE GOING TO BE ESTIMATED
alp = 0.75; % for C matrix
bet = 0.75; % for Q matrix


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
    for j=1:size(kin,2)

        % get current position and decoded velocity
        pos = kin(1:3,j);
        recorded_vel = TrialData.IntendedCursorState(:,j);
        %recorded_vel = kin(:,j);
        
        % if inside target, norm velocity = 0 else compute it
        Cursor.State = pos;
        Cursor.Center = TrialData.Params.Center;        
        out = InTargetRobot3D(Cursor,target_pos,...
            TrialData.Params.RobotTargetRadius);
        if out>0
            vel_mag= 0;
        else
            vel_mag = norm(recorded_vel(4:6));
        end
        
        % ReFit -> rotate decoded vel vector towards target
        new_vel_vector = target_pos' - pos;
        new_vel_vector = vel_mag*((new_vel_vector)/norm(new_vel_vector));
        xint = [xint new_vel_vector];
    end
    
    % delta, beta and hG
    fidx = [129:256 513:640 769:896]; 
    neural_features = neural_features(fidx,:);    

    % pooling
    if ch_pooling == true
        neural_features = pool_features(neural_features,TrialData.Params.ChMap);  
    end

    % store
    neural = [neural neural_features];
    kinematics = [kinematics xint];     
end


% STATES
X=kinematics;
X(end+1,:)=1;

% MEASUREMENTS
Y=neural;

% A matrix
A = KF_robot.A;

% W matrix
W = KF_robot.W;

% least squares estimates of velocity C
C = Y/X;
C = [zeros(size(C,1),3) C];

% get estimates of Q 
q = Y-C(:,4:end)*X;
Q= cov(q');
%Q = (1/size(Y,2)-1)* (q*q');


% SMOOTH BATCH UPDATE OF C AND Q
C = alp*KF_robot.C + (1-alp)*C;
Q = bet*KF_robot.Q + (1-bet)*Q;

% compute Kalman Gain and Error Covariance 
P = rand(7);
K = randn(size(C'));
chk=1;chk1=1;
counter=0;
norm_val=[];
norm_val1=[];
iter=1;
while chk>1e-10 || chk1>1e-10 
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
    if iter==1000
        break
    end
end

% store results
KF_robot.C = C;
KF_robot.Q = Q;
KF_robot.P = P;
KF_robot.K = K;

save(fullfile('clicker','kf_robot_params'),'KF_robot','-v7.3','-nocompression');
disp('Smooth batch of ReFit Kalman Filter Done')



