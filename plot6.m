dates = {'20210312'};

for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = ['/media/sarah/VICTOR-2/bravo1/', date, '/GangulyServer/', date, '/Robot/'];
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
   for fileInd = 1:numel(files)
    fn = [blkDir, '/', files{fileInd}];
% 
       load(fn)
       if TrialData.TargetID < 7
       unpackTD;
       %numTrials(TD.TargetID) = numTrials(TD.TargetID) + 1;
      figure(3)
      subplot(1,2,1)
        hold on
        axis equal
%            plot3(pos(1,:)', -pos(2,:)', pos(3, :)')
       plot3(pos(1,end)', -pos(2,end)', pos(3, end)',  'ko', 'markersize', 5, 'markerfacecolor', c)
       plot3(pos(1,:)', -pos(2,:)', pos(3, :)', 'color', c, 'linewidth', 2)
%             plot3(pos(1,:)', -pos(2,:)', pos(3, :)', 'linewidth', 1)
       plot3(0,0,0, 'ko', 'markerfacecolor', 'g', 'markersize', 5)
       plot3(TD.TargetPosition(1), -TD.TargetPosition(2), TD.TargetPosition(3), 'o', 'color', 'k', 'markersize', 12, 'markerfacecolor', c)
       a = 375;
       xlim([-a, a])
       ylim([-a, a])
       zlim([-a, a])
           view([30,30 ])
       grid on
       end
   end
end
end
%%
figure(3)
       a = 350;
       xlim([-a, a])
       ylim([-a, a])
       zlim([-a, a])
       
       
       subplot(1,2,1)
       title('20210312')
              
       subplot(1,2,2)
       title('20210319')