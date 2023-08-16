
clear all

SaveMovie = 0;
    rate = 0.2;
basedir     = '/media/sarah/Data/bravo1/';
dates       = {'20230811'};

load('/media/sarah/Data/bravo1/20230811/GangulyServer/20230811/MultiTargets2D/105924/BCI_Fixed/Data0001.mat')
Params = TrialData.Params;
%% PTB

Screen('Preference', 'SkipSyncTests', 1)
% [Params.WPTR, ScreenRectangle] = Screen('OpenWindow', max(Screen('Screens')), 0);
[Params.WPTR, ScreenRectangle] = Screen('OpenWindow', 0, 0, [0 0 1500 1200])
Screen('TextFont',Params.WPTR, 'Arial');
Screen('TextSize',Params.WPTR, 28);




%%

arrow = [1,0; 0, -1; -1,0; 0, 1];



%%
for dateInd = 1:numel(dates)
    date = dates{dateInd};
    taskDir = [basedir, date, '/GangulyServer/', date, '/MultiTargets2D/'];
    tmp = dir(taskDir);
    blockDirs = {tmp.name};
    blockDirs = blockDirs(3:end);

for blockInd = 2
    
    blkDir = [taskDir, blockDirs{blockInd}, '/BCI_Fixed'];
    tmp = dir(blkDir);
    files = {tmp.name};
    files = files(3:end);
    
    if  SaveMovie
        moviePtr = Screen('CreateMovie', Params.WPTR, blockDirs{blockInd},[], [],  Params.ScreenRefreshRate);

    end


    for fileInd = 1:length(files)

        
    fn = [blkDir, '/', files{fileInd}];
    load(fn)
    TD = TrialData;

    Params = TD.Params;

    Params.Center = [mean(ScreenRectangle([1,3])),mean(ScreenRectangle([2,4]))+200];


%% Fix CursorState

find(TD.TaskState < 3 )
st = find(TD.TaskState == 3,1);
en = find(TD.TaskState == 4,1)-1;
len  =length(find(TD.TaskState == 3));
ind = st:2:(st+(len-1)*2);

ind_tot = [find(TD.TaskState < 3 ), ind, ind(end)+2:length(TD.CursorState)];

TD.CursorState = TD.CursorState(:,ind_tot);

%% Set up trial
ReachTargetPos  = TD.TargetPosition;
TargetID        = TD.TargetID; % Target that cursor is in, 0 for no targets
TargetNum       = TD.Params.TargetList(TargetID,1);
StartPosNum     = TD.Params.TargetList(TargetID,2);

% Output to Command Line
fprintf('\nTrial: %i\n',TD.Trial)

%% Instructed Delay

bin_inds = 1:find(TD.TaskState == 2,1)-1;
for bin = bin_inds

    if SaveMovie
        Screen('AddFrameToMovie', Params.WPTR)
    end
       
            figure(1) 
        clf
        subplot(2,1,1)
        bar(TD.Belief(:,bin))
        hold on
        ylim([0, 1.0])
        title("Belief")
       subplot(2,1,2)
        hold on
        axis equal
        xlim([-1, 1])
        ylim([-1,1])
        dir = TD.FilteredClickerState(bin);
        if dir == 0
            plot(0,0, 'ko', 'markerfacecolor', 'k', 'markersize', 16)
        elseif dir < 5
            l = [0,0; arrow(dir,:)];
            plot(l(:,1), l(:,2),'r', 'linewidth', 2)
        elseif dir == 5
            plot(0,0, 'mo', 'markerfacecolor', 'm', 'markersize', 16)
        end
        title("Decode")
        pause(rate)   
end % only complete if no errors

%% Cue time

    % display targets
    for k = 1:Params.NumTargets
     TargetPos = Params.ReachTargets(k,:);
    
     Rect = Params.TargetRect; % centered at (0,0)
     Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
     Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos

    Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
     TargetPos = [0, -200];
    
    end
    
     TargetPos = Params.ReachTargets(TargetNum,:);
     Rect = Params.TargetRect; % centered at (0,0)
     Rect([1,3]) =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
     Rect([2,4]) =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
    Screen('FillOval', Params.WPTR, [255,0,0],  Rect)
    
    
    Screen('DrawingFinished', Params.WPTR);
    Screen('Flip', Params.WPTR);

    if SaveMovie
        Screen('AddFrameToMovie', Params.WPTR)
    end
    

bin_inds = find(TD.TaskState == 2,1):find(TD.TaskState == 3,1)-1;
for bin = bin_inds
      figure(1) 
        clf
        subplot(2,1,1)
        bar(TD.Belief(:,bin))
        hold on
        ylim([0, 1.0])
        title("Belief")
       subplot(2,1,2)
        hold on
        axis equal
        xlim([-1, 1])
        ylim([-1,1])
        dir = TD.FilteredClickerState(bin);
        if dir == 0
            plot(0,0, 'ko', 'markerfacecolor', 'k', 'markersize', 16)
        elseif dir < 5
            l = [0,0; arrow(dir,:)];
            plot(l(:,1), l(:,2),'r', 'linewidth', 2)
        elseif dir == 5
            plot(0,0, 'mo', 'markerfacecolor', 'm', 'markersize', 16)
        end
        title("Decode")
        pause(rate)   
       
end % only complete if no errors
        
%% Go to reach target

bin_inds = find(TD.TaskState == 3,1):find(TD.TaskState == 4,1)-1;

for bin = bin_inds
    
    Cursor.State = TD.CursorState(:,bin);
    
    inTarget = InTargetMulti2D(Cursor,Params.ReachTargets);

    figure(1)
    

        % draw targets
        for k = 1:Params.NumTargets
         TargetPos = Params.ReachTargets(k,:);
         Rect           = Params.TargetRect; % centered at (0,0)
         Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [200,200,200],  Rect)
         TargetPos = [0, -200];
        end

         % draw trial target color
         TargetPos      = Params.ReachTargets(TargetNum,:);
         Rect           = Params.TargetRect; % centered at (0,0)
         Rect([1,3])    =  Rect([1,3]) +  TargetPos(1) + Params.Center(1); % add x-pos
         Rect([2,4])    =  Rect([2,4]) +  TargetPos(2) + Params.Center(2); % add y-pos
         Screen('FillOval', Params.WPTR, [255,0,0],  Rect)
          
        % color target if cursor inside
        if any(inTarget)
            StartTargetPos      = Params.ReachTargets(find(inTarget),:);
            StartRect           = Params.TargetRect; % centered at (0,0)
            StartRect([1,3])    = StartRect([1,3]) + StartTargetPos(1) + Params.Center(1); % add x-pos
            StartRect([2,4])    = StartRect([2,4]) + StartTargetPos(2) + Params.Center(2); % add y-pos

            Screen('FillOval', Params.WPTR, [0,200,0], StartRect)
        end
                 
        CursorRect                  = Params.CursorRect;
        CursorRect([1,3])           = CursorRect([1,3]) + Cursor.State(1) + Params.Center(1); % add x-pos
        CursorRect([2,4])           = CursorRect([2,4]) + Cursor.State(2) + Params.Center(2); % add y-pos
        
        
        assist = norm(TD.AssistVel(:,bin))>0;
        % save               
        if assist && Params.ChangeAssistColor && Params.Assist
            Screen('FillOval', Params.WPTR, [255,0,255], CursorRect)
        else
            Screen('FillOval', Params.WPTR, [0,0,255], CursorRect)
        end

        Screen('DrawingFinished', Params.WPTR);
        Screen('Flip', Params.WPTR);
        if SaveMovie
            Screen('AddFrameToMovie', Params.WPTR)
        end

           figure(1) 
        clf
        subplot(2,1,1)
        bar(TD.Belief(:,bin))
        hold on
        ylim([0, 1.0])
        title("Belief")
       subplot(2,1,2)
        hold on
        axis equal
        xlim([-1, 1])
        ylim([-1,1])
        dir = TD.FilteredClickerState(bin);
        if dir == 0
            plot(0,0, 'ko', 'markerfacecolor', 'k', 'markersize', 16)
        elseif dir < 5
            l = [0,0; arrow(dir,:)];
            plot(l(:,1), l(:,2),'r', 'linewidth', 2)
        elseif dir == 5
            plot(0,0, 'mo', 'markerfacecolor', 'm', 'markersize', 16)
        end
        title("Decode")
        
        pause(rate)   
    
end
%% Inter trial interval
% blank screen at end of trial but continue collecting data

Screen('DrawingFinished', Params.WPTR);
Screen('Flip', Params.WPTR);

bin_inds = 1:find(TD.TaskState == 2,1)-1;
for bin = bin_inds   
       
    if SaveMovie
        Screen('AddFrameToMovie', Params.WPTR)
    end
               figure(1) 
        clf
        subplot(2,1,1)
        bar(TD.Belief(:,bin))
        hold on
        ylim([0, 1.0])
        title("Belief")
       subplot(2,1,2)
        hold on
        axis equal
        xlim([-1, 1])
        ylim([-1,1])
        dir = TD.FilteredClickerState(bin);
        if dir == 0
            plot(0,0, 'ko', 'markerfacecolor', 'k', 'markersize', 16)
        elseif dir < 5
            l = [0,0; arrow(dir,:)];
            plot(l(:,1), l(:,2),'r', 'linewidth', 2)
        elseif dir == 5
            plot(0,0, 'mo', 'markerfacecolor', 'm', 'markersize', 16)
        end
        title("Decode")

end
%%
    end


if SaveMovie
    Screen('FinalizeMovie', moviePtr);
end
end
end