function [p] = doGaborSteeringExample(subj, stimMode, viewdist)
% function [p] = yourPackage.doExperimentTemplate(subj, stimMode, viewdist)
% 
% Example for modular PLDAPS experimental design
% - based on modularDemo.doRfPos_gabGrid
% - "do" in filename is convention for PLDAPS experiment wrapper file that one would call
%   from the command window to run an experiment
%   e.g.  >> p = modularDemo.doRfPos_gabGrid("me", "gabors", 70);
% - 
% 
% INPUTS:
%   [subj]      Subject/session identifier  ("string", default: "test")
%   [stimMode]  Stimulus type               ("string", default: "gabors")   (...no other modes currently defined)
%   [viewdist]  
% 
% 
%   [pss] is the PLDAPS settings struct that is used to initialize modules
%   & parameters for this experimental session, and makeup the standard [p.trial]
%   structure that is integral to PLDAPS.
% 
% ----------------
% USAGE:
%   No inputs necessary to run demo from command window:
%   >> p = modularDemo.doRfPos_gabGrid
% 
% 
%   TRACKING CALIBRATION
%   To adjust eye/mouse tracking calibration:
%   - pause the experiment by pressing [p] key during a trial
%   - From the command window, start a tracking calibration trial:
%     >> pds.tracking.runCalibrationTrial(p)
%   - Follow directions printed in command window; briefly:
%     - [0] to reveal first fixation point
%     - [spacebar] to record fixation & advance to next point
%     - [u] to update calibration once sufficient points (>=10) have been recorded
%     - [p] to exit calibration, save to file, & return to pause state
%   
% 
% 2018-xx-xx TBC  Wrote it for RF mapping
% 2020-10-xx TBC  Updated to use consistent OpenGL rendering coordinates at screen center
% 2020-11-05 TBC  Cleaned & commented for modular tutorial
% 


KbName('UnifyKeyNames');

if nargin<1 || isempty(subj)
    subj = "test"; % "string" input preferred over 'char'
end
% stimMode == statename of stimulus module of interest
if nargin<2 || isempty(stimMode)
    stimMode = 'gabors';
end

% Set this module as the active stimulus module
%  (tells pmBase where to look for frame limits...incomplete implementation)
pss.pldaps.modNames.currentStim = {stimMode};

    
if nargin<3 || isempty(viewdist)
    viewdist = 57.29; % 57.29 == (1cm == 1deg)
end


pss.pldaps.pause.preExperiment = 0;

% play sound as/with reward feedback
pss.sound.use = 1;


% TUTORIAL TIP:
%   When learning/debugging PLDAPS, its helpful to set debug points inside your
%   modules so that you can examine how elements of your experiment operate/interact.
%   Overriding the default pldaps.trialMasterFunction ('runModularTrial') with the
%   following line:
% pss.pldaps.trialMasterFunction = 'runModularTrial_frameLock';
%   will allow you to manually step through your code without trial time elapsing.
%   -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
%   -----DO NOT use the  _frameLock  variant in your normal experiments!--------------
%   Since this _frameLock version completely breaks PLDAPS time keeping accuracy, a prominent
%   warning will be displayed in the command window when this trial function is used.
%   When running an proper experiment, best practice is to let .trialMasterFunction
%   inherit the default value from pldapsClassDefault.m by not setting it at all.
%   -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
% 

pss.newEraSyringePump.use = false;
pss.newEraSyringePump.refillVol = 40;
pss.behavior.reward.defaultAmount = 0.16;

fixRadius = 1; %[deg], applied to fixation module below


%% Eyepos & Eyelink
pss.eyelink.use = false;

pss.eyelink.useAsEyepos = pss.eyelink.use;
% make .mouse param respect eyelink use state
pss.mouse.useAsEyepos = ~pss.eyelink.useAsEyepos;
pss.pldaps.draw.eyepos.use = true;


%% Module inventory:
% -100: pldaps default trial function
%  1:   fixation
%  2:   base timing module
%  10:  gabor stim drawing module


%% display settings
pss.display.viewdist = viewdist;  % 57.29 == (1cm == 1deg)
pss.display.ipd = 6.5;  % human == 6.5;  macaque == 3.25;
pss.display.useOverlay = 1;

% pss.display.screenSize = [];
% pss.display.scrnNum = 0;

pss.display.stereoMode = 0; % 0==mono 2D; 3-or-4==freeFuse 3D;  See Screen('OpenWindow?')

pss.display.useGL = 1;
pss.display.multisample = 2;


%% (-100) pldaps default trial function
sn = 'pdTrialFxn'; % modName
pss.(sn) = pldapsModule('modName',sn, 'name','pldapsDefaultTrial', 'order',-100);


%% (1) fixation module
% see help  modularDemo.pmFixDot
sn = 'fix';
pss.(sn) =  pldapsModule('modName',sn, 'name','modularDemo.pmFixDot', 'order',1);

pss.(sn).use = true;
pss.(sn).on = true;
pss.(sn).mode = 0;          % eye position limit mode (0==pass/none, 1==square, 2==circle[euclidean]);  mode==2 strongly recommended;
pss.(sn).fixPos = [0 0];    % fixation xy in vis.deg, z in cm;    %NOTE: Recommended to leave z-position empty to ensure fixation is always rendered at current .display.viewdist
pss.(sn).fixLim = fixRadius*[1 1]; % fixation window limits (x,y)  % (visual degrees; if isscalar, y==x; if mode==2, radius limit; if mode==1, box half-width limit;)
pss.(sn).dotType = 12;      % 2 = classic PTB anti-aliased dot, 3:9 = geodesic sphere (Linux-only), 10:22 3D sphere (mercator; slow) (see help text from pmFixDot module for info on extended dotTypes available)
if pss.(sn).dotType <=2
    % pixels if classic PTB fixation dot
    pss.(sn).dotSz = 5; 
else
    % visual degrees of OpenGL dot
    pss.(sn).dotSz = 0.5; % vis. deg
end

% set this module as the active fixation module
% -- This is used to get/assign/update current .eyeX, .eyeY, .deltaXY positions
pss.pldaps.modNames.currentFix = {sn};


%% (2) base trial timing & behavioral control module
sn = 'pmBase';
pss.(sn) =  pldapsModule('modName',sn, 'name','templates.pmTemplate_track2fix', 'order',2);

stimDur = 3.6;
pss.(sn).stateDur = [NaN, 0.24, stimDur, NaN];

pss.(sn).waitForGo = 2;     % wait for button press (right shift key) before beginning trial stimulus presentation
pss.(sn).waitForResp = 1; % wait for (.hasResponse == true) || (.stateDur time expired)


%% (10) drifting gabor module:  glDraw.pmMatrixGabs.m

    %% Initialization
    % Initialize matrixModule(s) inside your "do" experiment wrapper by first setting up
    % [sn] as a normal struct field of your PLDAPS settings struct [pss]:

    %% Setup shared matrixModule parameters & components
    sn = 'gabors';

    % Base params field (non-module) stores params shared across all matrixModule instances
    pss.(sn).stimCtr = [0, 0, 0]; % stimulus center location; visual degrees
    % - The .stimPos of stimuli in each matrixModule will be relative to this common center location
    pss.(sn).gabTf = 3; % gabor temporal frequency
    pss.(sn).gabFwhm = 2; % full-width half-max size of each gabor (visual degrees)

    pss.(sn).centerOnScreen = true; % flag to use procedural gabor creation that has it's origin at screen center (rather than upper left corner)
    pss.(sn).trackMouse = 1;        
    pss.(sn).stimOffset = [0 0];

    % Create a temporary module that will be used to create each of your matrixModules
    tmpModule = pldapsModule('modName',sn, 'name','templates.pmTemplate_matrixModuleGabor', 'matrixModule',true, 'order',10); 

    % Apply any adjustments or defaults necessary in the matrixModule
    tmpModule.use = true;   % whether or not PLDAPS state machine should use the module when iterating through experiment states
                            % - important during basic initializaiton states like .experimentPreOpenScreen, 
    tmpModule.on = false;   % whether or not the module should be seen by your behavioral state machine (see modularDemo.pmBase.m)
  
    % Populate shared stimulus parameters 
    tmpModule.gabContrast = 1;
    tmpModule.ngabors = 1;
    % Some parameters are direclty copied over from the stimulus struct:
    % - NOTE that its just .(sn) during setup, but within the module code this shared location will become .(snBase)
    tmpModule.centerOnScreen = pss.(sn).centerOnScreen;
    tmpModule.gabTf = pss.(sn).gabTf; % drift rate (Hz)
    tmpModule.gabFwhm = pss.(sn).gabFwhm;
    tmpModule.gabSf = 1.3/tmpModule.gabFwhm;    % V1 median bandwidth= ~1.4 (== 1.0/fwhm) --De Valois 1982

    % Initialize parameters to be updated by condMatrix prior to each trial
    tmpModule.pos = zeros(2, tmpModule.ngabors); % [x,y] position of stimulus; visual degrees **relative to .stimCenter** (see code below)
    tmpModule.dir = zeros(1, tmpModule.ngabors);
    
    % %     % Select module stimulus type
    % %     % - Not necessary, but sometimes useful to have multiple 'flavors' of a particular module
    % %     % - Makes experiment code more robust/universal
    % %     % - Biggest benefit is ability to use common analysis code across experiments with the same general stimulus elements (i.e. gabors, dots, etc)
    % %     tmpModule.type = 'steerX';    %'polarTrack';
    
    % Stimulus onset timing & n-reps per trial
    ncopies = 1;
    tmpModule.isi = .0; % blank interstimulus interval between each matrix module
    % stimDur==total "stimulus on" duration within the trial
    % - should be setup while initializing your behavioral module (e.g. modularDemo.pmBase.m)
    stimModuleDur = stimDur/ncopies;


    %% Create matrixModules
    % Make duplicate/indexed stim modules for each repetition w/in a trial
    matrixModNames = {};
    for i = 1:ncopies
        mN = sprintf('%s%0.2d', sn, i); % indexed module name
        pss.(mN) = tmpModule;   % apply the temporary module
        pss.(mN).stateFunction.modName = mN; % Update it's internalized name field
        
        % timing: each module [onset, offset] time in sec (relative to STIMULUS state start)
        % - timing imposed by behavioral control module (e.g. modularDemo.pmBase)
        basedur = (i-1)*stimModuleDur;
        pss.(mN).modOnDur = [0, stimModuleDur-tmpModule.isi] +basedur; % subtract isi from offset

        % Compile module names [for use later within experiment wrapper code (...if needed)]s
        matrixModNames{i} = mN;
    end

    
%% Create PLDAPS object for this experimental session
p = pldaps(subj, pss);

% 
% LEGACY NOTE:
%   use of pdsDefaultTrialStructure.m is no longer recommended.
%       -- TBC 2020
% 


%% Define parameters of condition matrix

do3d = p.trial.display.stereoMode>0;

% Drift directions
dinc = 90; % direction increment
dirs = dinc:dinc:360; % motion directions

% starting stimulus location
xs = [-8, 8];
ys = [-4, 4];

% make a fully crossed matrix
[xx, yy, dd] = ndgrid( xs, ys, dirs);


% Create cell struct of module fields to be defined by each condition
% name of pldaps module that contains experimental manipulation
% ...how could we better know/extract this when needed?
% sn = stimMode;

% Make the condMatrix!
c = cell(size(xx)); % *** maintain same shape as condition matrix source values

for i = 1:numel(xx)
    % Set up conditions
    c{i}.stimPos    = [ xx(i), yy(i)]; % stim position relative to center
    c{i}.dir    = dd(i) * ones(tmpModule.ngabors, 1+do3d); % separate dir for left & right eye
    
end


%% Generate condMatrix & add to PLDAPS
p.condMatrix.conditions = c;
% 
% 'randMode' & 'nPasses' are the primary name-value pairs used to control the condMatrix:
% 
%   [nPasses]   (default: inf)
%               # of full passes through condition matrix to be completed 
%               - inf will run until experiment is manually quit ([q] key pressed)
%               - presentationw will continue until nPasses are COMPLETE,
%                 padding out the final trial with additional stimuli [from the next pass]
%                 (e.g. a 4-by-4 condition matrix with 6 stimuli per trial,
%                  nPasses=1 will comprise 4 trials, with a repeat of the first two matrix
%                  conditions; filling out the last trial.
%                 
% 
%   [randMode]  (default: 0)
%               Randomize order of upcoming pass through condition matrix
%               --Simple--
%               0 == no randomization; walk through matrix in sequential [column-major] order
%               1 == randomize across all dimensions
%               2 == randomize within columns w/ Shuffle.m
%               3 == randomize within rows (**primarily only for 2D condMatrix, 3D hacky, >3D will error)
%               --Indexed--
%               Index dimensions to be randomized:
%                 - positive randMode values will randomize each element of dimension
%                 - zero randMode values will do nothing
%                 - negative randMode values will shuffle the order of that dimension,
%                   while maintaining other dimension(s)
%                 Ex: condMatrix of [xpos, ypos, direction]
%                     - randMode [1,2,3] will have same effect as randMode [1]
%                     - randMode [1,2,-3] will present one random direction at each x,y position randomly,
%                       before advancing to a new random direction, etc, etc
%                     - randMode [3,0] will present x,y positions in order, with a different random
%                       direction at each sequential location
%                       (* here the 0 distinguishes it from simple randMode [3])
%                       

p.condMatrix = condMatrix(p, 'randMode',1, 'nPasses',inf);
% 
% p.condMatrix = condMatrix(p, 'randMode',[2], 'nPasses',2);

% p.condMatrix = condMatrix(p, 'randMode',[3,0], 'nPasses',1);

% p.condMatrix = condMatrix(p, 'randMode',[0], 'nPasses',1);

%
% LEGACY NOTE: 
%   !! Leave p.conditions empty !! everything is in .condMatrix now
%   old-style conditions definitions should still work through PLDAPS v4.x, but
%   testing/support for that method is waining & will eventually error.
%       --TBC 2020
% 


%% Run it!!
p.run;

% % Create a summary plot of display timing accuracy
% pds.plotTiming(p);
%     % 
%     % **NOTE: This would NOT be especially meaningful if using
%     % the '_frameLock' mode described at top of this demo, but is
%     % a good starting point for assessing accuracy during standard
%     % PLDAPS usage.
%     % (i.e. w/default:  .pldaps.trialMasterFunction = 'runModularTrial';)
%     % 

end
