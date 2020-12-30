function p = pdsDefaultTrialStructure(p, stimulusname)
% p = pdsDefaultTrialStructure(p)
% pdsDefaultTrialStructure sets up the p struct for displaying all the
% standard task parameters used in huk lab
% 
% UPDATED GUIDANCE (circa 2020)
% -----------------------
% Use of this setup function is no longer necessary/recommended.
% Old code excessively relied on a hardcoded  "p.trial.stimulus"  element
% of the pldaps object.
% 
% Modern coding approaches in PLDAPS lean on a more flexible system
% of "pldaps modules" which expect the [p.trial.(sn)] name to be passed
% as an input string to the module code, e.g.
%   function p = modularDemo.pmBase(p, state, sn)
%   % ...
%   %   switch state
%   %       case p.trial.pldaps.trialStates.frameUpdate
%   %   ...
% 
% For an example of a modular PLDAPS replacement for much of what this
% function did inside of the [p.trial.stimulus] field,
% see also: modularDemo.pmBase
% 
% ---
% 2020-12-30  TBC  Added comments clarifying proper replacement
% 


% Discourage use  --TBC 2020
fprintf(2, fprintLineBreak('~'))
fprintf(2, '~!~\t\tWARNING  %s.m  WARNING  \t\t~!~\n', mfilename);
fprintf(2, '~!~\t\t   Usage of this function is deprecated \t\t\t~!~\n');
fprintf('Please update your code for use with modular PLDAPS style\n')
fprintf(2, fprintLineBreak('~'))


if nargin<2
    stimulusname='stimulus';
end

% 12/2013 jly   Wrote it
p.defaultParameters.good = 1;

% p.defaultParameters.(stimulusname).photodiode.use = 1;
% p.defaultParameters.(stimulusname).photodiode.location = 2;
% p.defaultParameters.(stimulusname).photodiode.frames = 10;


if ~isfield(p.defaultParameters,'pldaps.finish')
    p.defaultParameters.pldaps.finish = inf;
end

p.defaultParameters.(stimulusname).randomNumberGenerater = 'mt19937ar';


% Setup Timings
%-------------------------------------------------------------------------%
% Reward time: is time that solenoid is opened for. set to 100 miliseconds
% for mapping trials.
p.defaultParameters.(stimulusname).rewardTime = .1;
p.defaultParameters.(stimulusname).rewardWait = 0;
p.defaultParameters.(stimulusname).breakFixPenalty = 2;
p.defaultParameters.(stimulusname).jitterSize = .5;
% fixation
p.defaultParameters.(stimulusname).preTrial     = .5;
p.defaultParameters.(stimulusname).fixWait      = 4;
p.defaultParameters.(stimulusname).fixHold      = 1;

% targets
p.defaultParameters.(stimulusname).targWait   = 1.5;
p.defaultParameters.(stimulusname).targHold   = 0.5;
p.defaultParameters.(stimulusname).targOnset  = [0.1 0.1];
p.defaultParameters.(stimulusname).targDuration = [2 .2];

% Colors
%-------------------------------------------------------------------------%
p = defaultColors(p);

% Bits
%-------------------------------------------------------------------------%
p = defaultBitNames(p);

% dot sizes for drawing
p.defaultParameters.(stimulusname).eyeW      = 8;    % eye indicator width in pixels
p.defaultParameters.(stimulusname).fixdotW   = 8;    % width of the fixation dot
p.defaultParameters.(stimulusname).targdotW  = 8;    % width of the target dot
p.defaultParameters.(stimulusname).cursorW   = 8;   % cursor width in pixels

% States
%-------------------------------------------------------------------------%
p.defaultParameters.(stimulusname).states.START     = 1;
p.defaultParameters.(stimulusname).states.FPON      = 2;
p.defaultParameters.(stimulusname).states.FPHOLD    = 3;
p.defaultParameters.(stimulusname).states.CHOOSETARG = 4;
p.defaultParameters.(stimulusname).states.HOLDTARG     = 5;
p.defaultParameters.(stimulusname).states.BREAKFIX  = 7;
p.defaultParameters.(stimulusname).states.TRIALCOMPLETE = 6;



