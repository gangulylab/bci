function KF = UpdateCursor(Params,Neuro,TaskFlag,TargetPos,KF)
% KF = UpdateCursor(Params,Neuro)
% Updates the state of the cursor using the method in Params.ControlMode
%   1 - position control
%   2 - velocity control
%   3 - kalman filter position & velocity
% 	4 - kalman filter velocity
%
% Cursor - global structure with state of cursor [px,py,vx,vy,1]
% TaskFlag - 0-imagined mvmts, 1-clda, 2-fixed decoder
% TargetPos - x- and y- coordinates of target position. used to assist
%   cursor to target
% KF - kalman filter struct containing matrices A,W,P,C,Q

global Cursor

% query optimal control policy
Vopt = OptimalCursorUpdate(Params,TargetPos);

if TaskFlag==1, % do nothing during imagined movements
    return;
end

% find vx and vy using control scheme
switch Cursor.ControlMode,
    case 1, % Move to Mouse
        X0 = Cursor.State;
        [x,y] = GetMouse();
        vx = ((x-Params.Center(1)) - Cursor.State(1))*Params.UpdateRate;
        vy = ((y-Params.Center(2)) - Cursor.State(2))*Params.UpdateRate;
        
        % update cursor
        Cursor.State(1) = x - Params.Center(1);
        Cursor.State(2) = y - Params.Center(2);
        Cursor.State(3) = vx;
        Cursor.State(4) = vy;
        
        % Update Intended Cursor State
        X = Cursor.State;
        Vcom = (X(1:2) - X0(1:2))*Params.UpdateRate; % effective velocity command
        Cursor.IntendedState = Cursor.State; % current true position
        Cursor.IntendedState(3:4) = Vopt; % update vel w/ optimal vel
        
    case 2, % Use Mouse Position as a Velocity Input (Center-Joystick)
        X0 = Cursor.State;
        [x,y] = GetMouse();
        vx = Params.Gain * (x - Params.Center(1));
        vy = Params.Gain * (y - Params.Center(2));
        
        % assisted velocity
        if Cursor.Assistance > 0,
            Vcom = [vx;vy];
            Vass = Cursor.Assistance*Vopt + (1-Cursor.Assistance)*Vcom;
        else,
            Vass = [vx;vy];
        end
        
        % update cursor state
        Cursor.State(1) = Cursor.State(1) + Vass(1)/Params.UpdateRate;
        Cursor.State(2) = Cursor.State(2) + Vass(2)/Params.UpdateRate;
        Cursor.State(3) = Vass(1);
        Cursor.State(4) = Vass(2);
        
        % Update Intended Cursor State
        X = Cursor.State;
        Vcom = (X(1:2) - X0(1:2))*Params.UpdateRate; % effective velocity command
        Cursor.IntendedState = Cursor.State; % current true position
        Cursor.IntendedState(3:4) = Vopt; % update vel w/ optimal vel
        
    case {3,4}, % Kalman Filter Input
        X0 = Cursor.State; % initial state, useful for assistance
        %X0 =[ .1 .1 .1 .1 .1]';
        
        % Kalman Predict Step
        X = X0;
        if Neuro.DimRed.Flag,
            Y = Neuro.NeuralFactors;
        else,
            Y = Neuro.MaskedFeatures;
        end
        A = KF.A;
        W = KF.W;
        P = KF.P;
        C = KF.C;
        Q = KF.Q;
        
        % Kalman Predict Step
        X = A*X;
        P = A*P*A' + W;
        P(1:2,:) = zeros(2,5); % zero out pos and pos-vel terms
        P(:,1:2) = zeros(5,2); % innovation from refit
        
        % Kalman Update Step
        K = P*C'/(C*P*C' + Q);
        X = X + K*(Y - C*X);
        P = P - K*C*P;
        
        % Store Params
        Cursor.State = X;
        KF.P = P;
        KF.K = K;
        
        % assisted velocity
        Vcom = X(3:4); % effective velocity command
        if Cursor.Assistance > 0,
            % Vass w/ vector avg
            %Vass = Cursor.Assistance*Vopt + (1-Cursor.Assistance)*Vcom;
            
            if ~Params.DaggerAssist, % Vass w/ same speed
                norm_vcom = norm(Vcom);
                Vass = Cursor.Assistance*Vopt + (1-Cursor.Assistance)*Vcom;
                Vass = norm_vcom * Vass / norm(Vass);
            else, % Dagger Assist
                sample_optimal = rand(1)<Cursor.Assistance;
                if sample_optimal,
                    Vass = Vopt;
                else, % sample regular
                    Vass = Vcom;
                end
            end
            
            % update cursor state
            %Cursor.State(1) = X0(1) + Vass(1)/Params.UpdateRate;
            %Cursor.State(2) = X0(2) + Vass(2)/Params.UpdateRate;
            Cursor.State(3) = Vass(1);
            Cursor.State(4) = Vass(2);
        end
        
        % Update Intended Cursor State
        Cursor.IntendedState = Cursor.State; % current true position
        Cursor.IntendedState(3:4) = Vopt; % update vel w/ optimal vel
        
        % Apply Velocity Transform
        if Params.VelocityTransformFlag,
            [Vx,Vy] = VelocityTransform(Cursor.State(3),Cursor.State(4),Params.Gain);
            Cursor.State(3) = Vx;
            Cursor.State(4) = Vy;
        elseif Params.MaxVelocityFlag,
            speed = norm(Cursor.State(3:4));
            new_speed = max([speed,Params.MaxVelocity]);
            Cursor.State(3) = Cursor.State(3) * new_speed / speed;
            Cursor.State(4) = Cursor.State(4) * new_speed / speed;
        end
        
        % Update KF Params (RML & Adaptation Block)
        if KF.CLDA.Type==3 && TaskFlag==2,
            KF = UpdateRmlKF(KF,Cursor.IntendedState,Y,Params,TaskFlag);
        elseif KF.CLDA.Type==3 && TaskFlag==3 && Params.CLDA.FixedRmlFlag, % (RML & Fixed)
            KF = UpdateRmlKF(KF,Cursor.State,Y,Params,TaskFlag);
        end
        
end

% update effective velocity command for screen output
try,
    Cursor.Vcommand = Vcom;
catch,
    Cursor.Vcommand = [0,0];
end

% bound cursor position to size of screen
pos = Cursor.State(1:2)' + Params.Center;
pos(1) = max([pos(1),Params.ScreenRectangle(1)+10]); % x-left
pos(1) = min([pos(1),Params.ScreenRectangle(3)-10]); % x-right
pos(2) = max([pos(2),Params.ScreenRectangle(2)+10]); % y-left
pos(2) = min([pos(2),Params.ScreenRectangle(4)-10]); % y-right
Cursor.State(1) = pos(1) - Params.Center(1);
Cursor.State(2) = pos(2) - Params.Center(2);

end % UpdateCursor