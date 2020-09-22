function [ Params ] = SetKeyboardParams(Params, varargin)
% [ Params ] = SetKeyboardParams(Params)

p = inputParser;
p.addOptional('NoDisplay', false, @islogical);
parse(p, varargin{:});

KP = struct();

if p.Results.NoDisplay
    fprintf('No Display\n')
    Params = struct('Center', [0, 0], 'ReachTargetPositions', zeros(8, 2));
end

% Params you may want to change during an experiment
KP.Text.CorpusSize = 900;
KP.MaxUndo = 3;
KP.N_Back = 0;
% debugging
KP.Verbose = false;
% Targets
KP.TargetWidth = Params.OuterCircleRadius;
KP.TargetHeight = Params.OuterCircleRadius;
% TrialData.TargetPosition normally set in RunLoop.m
% KP.TargetPosition = Params.ReachTargetPositions + Params.Center; % old typing env
KP.TargetPosition = Params.ReachTargetPositions + Params.Center; % from new grid task setup
% Params.TargetRect normally set in GetParams.m
% NOTE: removing for grid task
KP.TargetRect = [-KP.TargetWidth / 2, -KP.TargetHeight / 2, KP.TargetWidth / 2, KP.TargetHeight / 2]; % Left, Top, Right, Bottom

KP.CharColor     = [255, 255, 0];
KP.WordColor     = [0, 150, 255];
KP.NextWordColor = [0, 255, 150];

%% Text
KP.Text.CharacterSets = {'A B','A B','A B','A B'};%,'A B','A B','A B','A B'};
KP.Text.CharDisplayOpts = {'FontSize', 25,...
                            'Offset', [10, 18],...
                             'Color', [170, 170, 170]};
KP.Text.WordDisplayOpts = {'FontSize', 32,...
                            'Offset', [10, 15],...
                             'Color', [170, 170, 170]};
KP.Text.WordBox.Title = 'Next Words';
KP.Text.WordBox.TextOpts = {'FontSize', 25,...
                            'Offset', [0, 0],...
                            'Color', [0, 0, 0],};
% KP.Text.NextWordBox.Title = 'Next Words';
% KP.Text.NextWordBox.TextOpts = {'FontSize', 25,...
%                             'Offset', [0, 0],...
%                             'Color', [0, 0, 0],};
KP.Text.WordSet = GetCorpus('N', KP.Text.CorpusSize);

%% Keyboard layout
% ix_F_Arrow = 1;
% ix_B_Arrow = 5;
ix_F_Arrow = [];
ix_B_Arrow = [];
KP.Pos = struct();
KP.Pos.ArrowTargets = KP.TargetPosition([ix_F_Arrow, ix_B_Arrow], :);
KP.Pos.ArrowLabels = {};
KP.Pos.TextTargets = KP.TargetPosition;
KP.Pos.F_Arrow = [];
KP.Pos.B_Arrow = [];
KP.Pos.TargetEdges = (repmat(KP.TargetPosition, 1, 2) + KP.TargetRect)';
[~, ix_top_targ] = min(KP.TargetPosition(:, 2));
[~, ix_right_targ] = max(KP.TargetPosition(:, 1));
KP.Pos.CharDisplay = KP.TargetPosition(ix_top_targ, :) - [0, KP.TargetHeight * 0.60];
KP.Pos.WordDisplay = KP.Pos.CharDisplay - [0, KP.TargetHeight * 0.2];

% End screen
ix_end_targets = [ix_F_Arrow, ix_B_Arrow, ix_top_targ];
KP.End.Targets = KP.TargetPosition(ix_end_targets, :);
KP.End.TargetLabels = {'CONTINUE', 'STOP', 'Reset'};
KP.End.TargetEdges = KP.Pos.TargetEdges(:, ix_end_targets);
KP.End.Color = [0, 255, 0;
                255, 0, 0;
                0, 150, 255]';

% Word Box
KP.ShowWordBox = false;
KP.Pos.WordBox.H = KP.TargetHeight * 1.5;
KP.Pos.WordBox.W = KP.TargetWidth;
KP.Pos.WordBox.O = [KP.TargetPosition(ix_right_targ, 1), KP.TargetPosition(ix_top_targ, 2)] + [KP.TargetWidth, 0];
KP.Pos.WordBox.Edges = [KP.Pos.WordBox.O, KP.Pos.WordBox.O + [KP.Pos.WordBox.W, KP.Pos.WordBox.H]];
% KP.Pos.WordBox.Color = KP.WordColor;
KP.Pos.WordBox.FirstEntry = KP.Pos.WordBox.O + [10, 10];
KP.Pos.WordBox.WordSpacing = [0, 40];

%% State
n_text = size(KP.Pos.TextTargets, 1);
n_arrow = size(KP.Pos.ArrowTargets, 1);
KP.State = struct();
KP.State.CurrentTargets = KP.TargetPosition;
KP.State.CurrentTargetEdges = KP.Pos.TargetEdges;
KP.State.InText = false(n_text, 1);
KP.State.InArrow = false(n_arrow, 1);
KP.State.InCurrent = false(size(KP.State.CurrentTargets, 1), 1);
KP.State.Mode = 'Character'; % or 'Word' - set during task
KP.State.Select = false; % indicates selection has been made
KP.State.SelectableText = KP.Text.CharacterSets;
KP.State.NText = n_text;
KP.State.NArrow = n_arrow;
KP.State.NextWordSet = KP.Text.WordSet(1:n_text); % NOTE: text here, not well organized
KP.State.CurrentColor = KP.CharColor;
KP.State.WordBoxColor = KP.WordColor;
KP.State.SelectedCharacters = {};
KP.State.SelectedWords = {};
KP.State.WordMatches = 1:length(KP.Text.WordSet);
KP.State.TargetID = 0;

% Override target color
Params.TargetsColor         = KP.State.CurrentColor; % all targets

% Keyboard State History
KP.History.State = {};
KP.History.InitState = KP.State;

Params.Keyboard = KP;
end  % function
