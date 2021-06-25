function [Neuro,KF,Params] = RunTask(Params,Neuro,TaskFlag,KF)
% Explains the task to the subject, and serves as a reminder for pausing
% and quitting the experiment (w/o killing matlab or something)

global Cursor 
Cursor.ControlMode = Params.ControlMode;

% Load Clicker
% f = load(fullfile('clicker',Params.DiscreteDecoder)); % load the parallel decoders
% Clicker.Model = f.model;
% Clicker.Func = @ (X) multistate_discrete(X,Clicker.Model,Params.MultiDecisionBoundary); % calling function
Clicker = []
switch TaskFlag
    case 1 % Imagined Movements

        Cursor.Assistance = Params.Assistance;
        Cursor.DeltaAssistance = 0;
        InstructedDelayTime = Params.InstructedDelayTime;
%         Params.InstructedDelayTime = .6;
        mkdir(fullfile(Params.Datadir,'Imagined'));
        
        % output to screen
        fprintf('\n\nImagined Movements:\n')
        fprintf('  %i Blocks (%i Total Trials)\n',...
            Params.NumImaginedBlocks,...
            Params.NumImaginedBlocks*Params.NumTrialsPerBlock)
        fprintf('  Saving data to %s\n\n',fullfile(Params.Datadir,'Imagined'))
        
        Neuro.DimRed.Flag = false; % set to false for imagined mvmts
        [Neuro,~,Params,~] = RunLoop(Params,Neuro,TaskFlag,fullfile(Params.Datadir,'Imagined'),[],[]);
        Params.InstructedDelayTime = InstructedDelayTime;
        
    case 2, % Control Mode with Assist & CLDA
        switch Params.ControlMode,
            case 1, % Mouse Position
                Instructions = [...
                    '\n\nMouse Position Control\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
            case 2, % Mouse Velocity
                Instructions = [...
                    '\n\nMouse Velocity Control\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
            case {3,4}, % Kalman Filter Decoder
                Instructions = [...
                    '\n\nKalman Brain Control (Calibration Mode)\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
                
                % Fit Dimensionality Reduction Params & Kalman Filter 
                % based on imagined mvmts
                Neuro.DimRed.Flag = Params.DimRed.Flag; % reset for task
                if Params.DimRed.Flag && Params.DimRed.InitAdapt,
                    Neuro.DimRed.F = FitDimRed(...
                        fullfile(Params.Datadir,'Imagined'),Neuro.DimRed,Params);
                    KF = FitKF(Params,...
                        fullfile(Params.Datadir,'Imagined'),0,KF,[],Neuro.DimRed.F);
                else, % no dim reduction
                    KF = FitKF(Params,...
                        fullfile(Params.Datadir,'Imagined'),0,KF);
                end
                
                
        end
        
        InstructionScreen(Params,Instructions);
        Cursor.Assistance = Params.Assistance;
        Cursor.DeltaAssistance = Params.CLDA.DeltaAssistance;
        mkdir(fullfile(Params.Datadir,'BCI_CLDA'));
        
        % output to screen
        fprintf('\n\nAdaptive Control: (%s)\n', Params.CLDA.TypeStr)
        fprintf('  %i Blocks (%i Total Trials)\n',...
            Params.NumAdaptBlocks,...
            Params.NumAdaptBlocks*Params.NumTrialsPerBlock)
        fprintf('  Assistance: %.2f\n', Cursor.Assistance)
        fprintf('  Change in Assistance: %.2f\n', Cursor.DeltaAssistance)
        fprintf('  Saving data to %s\n\n',fullfile(Params.Datadir,'BCI_CLDA'))
        
        [Neuro,KF,Params,Clicker] = RunLoop(Params,Neuro,TaskFlag,fullfile(Params.Datadir,'BCI_CLDA'),KF,Clicker);
        
    case 3, % Control Mode without Assist and fixed
        switch Params.ControlMode,
            case 1, % Mouse Position
                Instructions = [...
                    '\n\nMouse Position Control\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
            case 2, % Mouse Velocity
                Instructions = [...
                    '\n\nMouse Velocity Control\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
            case {3,4}, % Kalman Filter Decoder
                Instructions = [...
                    '\n\nKalman Brain Control\n\n'...
                    '\nAt any time, you can press ''p'' to briefly pause the task.'...
                    '\n\nPress the ''Space Bar'' to begin!' ];
                
                % reFit Kalman Filter based on intended kinematics during
                % adaptive block
                Neuro.DimRed.Flag = Params.DimRed.Flag; % reset for task
                if Neuro.CLDA.Type==1,
                    if Params.DimRed.Flag,
                        KF = FitKF(Params,...
                            fullfile(Params.Datadir,'BCI_CLDA'),1,KF,[],Neuro.DimRed.F);
                    else,
                        KF = FitKF(Params,...
                            fullfile(Params.Datadir,'BCI_CLDA'),1,KF);
                    end
                elseif Params.NumAdaptBlocks==0 || Neuro.DimRed.InitFixed,
                    if Params.DimRed.Flag,
                        Neuro.DimRed.F = FitDimRed(...
                            fullfile(Params.Datadir,'Imagined'),Neuro.DimRed,Params);
                        KF = FitKF(Params,...
                            fullfile(Params.Datadir,'Imagined'),0,KF,[],Neuro.DimRed.F);
                    else,
                        KF = FitKF(Params,...
                            fullfile(Params.Datadir,'Imagined'),0,KF);
                    end
                elseif Neuro.CLDA.Type==3, % update Qinv
                    KF.Qinv = inv(KF.Q);
                    KF.Lambda = Params.CLDA.FixedLambda;
                end
        end
        
        Cursor.Assistance = 0;
        Cursor.DeltaAssistance = 0;
        mkdir(fullfile(Params.Datadir,'BCI_Fixed'));
        
        % output to screen
        fprintf('\n\nFixed Control:\n')
        fprintf('  %i Blocks (%i Total Trials)\n',...
            Params.NumFixedBlocks,...
            Params.NumFixedBlocks*Params.NumTrialsPerBlock)
        fprintf('  Saving data to %s\n\n',fullfile(Params.Datadir,'BCI_Fixed'))
        
        [Neuro,KF,Params,Clicker] = RunLoop(Params,Neuro,TaskFlag,fullfile(Params.Datadir,'BCI_Fixed'),KF,Clicker);
        
end

end % RunTask
