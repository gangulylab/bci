dates = {'20210312'};

for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = ['/media/sarah/VICTOR-2/bravo1/', date, '/GangulyServer/', date, '/Robot/'];
    taskDir = ['/home/sarah/Downloads/Robot/'];
    tmp = dir(taskDir);
    blockDirs = {tmp.name};
    blockDirs = blockDirs(3:end);
% for blockInd = 1:numel(blockDirs)
for blockInd = 2
    blkDir = [taskDir, blockDirs{blockInd}, '/BCI_Fixed'];
    tmp = dir(blkDir);
    files = {tmp.name};
    files = files(3:end);
    numTrials = 0;
% for fileInd = 1
   for fileInd = 1
    fn = [blkDir, '/', files{fileInd}];
    fileInd
    
    
   
% 
% fn = '/media/sarah/VICTOR-2/bravo1/20210319/GangulyServer/20210319/Robot3DArrow/110258/BCI_Fixed/Data0001.mat';



load(fn)

Params = TrialData.Params;

Params.udp = udp("127.0.0.1", 5006);
fopen(Params.udp)
fwrite(Params.udp, [0,1,0])                  % reset robot
fwrite(Params.udp, [0,2,Params.UpdateRate])  % set update rate
fwrite(Params.udp, [0,3,Params.RobotMode])   % set robot mode
fwrite(Params.udp, [0,4,Params.RobotDirectionLines])   % set debug lines
%%

% fn = '/media/sarah/VICTOR-2/bravo1/20210219/GangulyServer/20210219/RobotRR/150252/BCI_Fixed/Data0001.mat';
load(fn)

TD = TrialData;
figure(1)
hold on
axis equal
plot3(TD.CursorState(1,:), TD.CursorState(2,:), TD.CursorState(3,:))

pause()
   end

        
    end


end
