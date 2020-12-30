function[p] = doFrameTimingTest(subj, stimMode, viewdist)

KbName('UnifyKeyNames');


%% Parse inputs
if nargin<1 || isempty(subj)
    subj = 'test';
end
% tstMode == statename of stimulus module of interest
if nargin<2 || isempty(stimMode)
    stimMode = 'dotBall';
end
% Set this module as the active stimulus module
%  (tells pmBase where to look for frame limits...incomplete implementation)
pss.pldaps.modNames.currentStim = {stimMode};

    
if nargin<3 || isempty(viewdist)
    viewdist = 101;% (101cm == 20ppd)    57.29; % 57.29 == (1cm == 1deg)
elseif ~isscalar(viewdist)
    % parse multiple viewdists into blocks
    vdists = viewdist(:);
    % set initial active viewdist to maximum
    % - ensures sufficient texture support
    viewdist = max(viewdist);
end


%% Toggle generic defaults
pss.tracking.use = true;
pss.pldaps.pause.preExperiment = true;

pss.pldaps.draw.grid.use = 0;
pss.pldaps.draw.eyepos.use = 0;

pss.newEraSyringePump.use = false;
% pss.newEraSyringePump.refillVol = 40;
% pss.newEraSyringePump.allowNewDiameter = true;
% pss.behavior.reward.defaultAmount = 0.02;

fixRadius = 2;

%% Eyepos & Eyelink
pss.eyelink.use = false;
pss.eyelink.useAsEyepos = pss.eyelink.use;
% make .mouse param respect eyelink use state
pss.mouse.useAsEyepos = ~pss.eyelink.useAsEyepos;

% 
% 
% DEBUG:  Mouse simulation mode
pss.tracking.useRawData = true;
% 
% 
% 


%% Module inventory:
% -100: pldaps default trial function
% -0.5: grbl (auto disabling)
%  0:   fixation
%  0.5: base timing module
%  1:   gabor stim drawing module
%  100: movie capture module (*causes frame drops when enabled, not for actual expt use*) 


%% display settings
pss.display.viewdist = viewdist; % 13 (7.7 diopters) ,20 (5 diopters), 31 cm (~3.25) diopters, 67 cm (~1.5 diopters)
pss.display.ipd = 3.25;
pss.display.useOverlay = 1;

pss.display.forceLinearGamma = false;
% pss.display.stereoMode = 0;
% pss.datapixx.rb3d = 0;

pss.display.useGL = 1;
pss.display.multisample = 2;

% OpenGL view settings
pss.display.obsPos = [0 0 0 0];
pss.display.fixPos = [0 0 pss.display.viewdist];
pss.display.upVect = [0 1 0];

pss.display.zNear = 0.1*pss.display.viewdist;
pss.display.zFar = 5*max([pss.display.viewdist, pss.display.fixPos(3), 100]);


%% (-100) pldaps default trial function
sn = 'pdTrialFxn'; % modName
pss.(sn) = pldapsModule('modName',sn, 'name','pldapsDefaultTrial', 'order',-100);


%% (-0.5) grbl module:  controls viewdist screen position
sn = 'grbl';
pss.(sn) = pldapsModule('modName',sn, 'name','grbl.grbl', 'order',-0.5,...
    'requestedStates',{'trialSetup','trialPrepare','experimentPreOpenScreen','experimentPostOpenScreen','experimentCleanUp'});
pss.(sn).use = false;


%% (1) fixation module
sn = 'fix';
pss.(sn) =  pldapsModule('modName',sn, 'name','visBasics.pmFixDot', 'order',1);

pss.(sn).use = false;
pss.(sn).on = false;
pss.(sn).mode = 0;
pss.(sn).fixPos = pss.display.fixPos(1:2); % [0,0,pss.display.viewdist];%pss.display.fixPos;
pss.(sn).fixLim = fixRadius*[1 1]; % fixation window limits (x,y) % 2.5; % outer radius in visual degrees
pss.(sn).dotType = 12; % 2 = classic PTB anti-aliased dot
pss.(sn).dotSz = 0.5; % pixels if dotType==2

% set this module as the active fixation module
% -- This is used to get/assign/update current .eyeX, .eyeY, .deltaXY positions
pss.pldaps.modNames.currentFix = {sn};


%% (2) base trial timing module
sn = 'pmBase';
pss.(sn) =  pldapsModule('modName',sn, 'name','glDraw.pmBase', 'order',2);

pss.(sn).waitForGoSignal = false;

fixDur = 3.6;
pss.(sn).stateDur = [NaN, 0.24, fixDur, NaN];


%% (1) drifting gabor module:  visBasics.pmDemoGabs.m

gridN   =  5*[1 1];    %[5, 9]; % 7*[1 1];    %7 * [1 1];    % number of samples in xy dimension
ipsiColumn = [1 0];

switch stimMode
    case 'dotBall'
        
        sn = stimMode;
        % Base params field (non-module) stores params shared across all matrixModule instances
        
        if 1 % smaller grid
%             pss.(sn).stimCtr = [11, -10, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [25, 25]; %30 *[1 1];
            
%             pss.(sn).stimCtr = [6.25, -1.75, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = 22 *[1 1];

% %             pss.(sn).stimCtr = [9, -6, 0]; % [x, y, zOffsetFromFixation]
% %             pss.(sn).gridSz  = 16 *[1 1];

            pss.(sn).stimCtr = [7.5, -12, 0]; % [x, y, zOffsetFromFixation]
            pss.(sn).gridSz  = 18 *[1 1];

%             pss.(sn).stimCtr = [4.3, -4.3, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = 11 *[1 1];

        else % searching grid   
            % FULL Right hemifield w/ ipsillateral column
%             % fwhm==    4.25     @ 7x7
%             pss.(sn).stimCtr = [30, 0, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [70,55];  %35 *[1 1];
            
%             % fwhm==    4.0     @ 7x7
%             pss.(sn).stimCtr = [21, -5, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [35,35];  %35 *[1 1];

% %             % MID-SIZED Lower Right hemifield w/ ipsillateral column
%             ipsiColumn = [1 0];
%             pss.(sn).stimCtr = [14.5, -7.5, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [35,30];  %35 *[1 1];
            
            % MID-SIZED Lower Right hemifield, slight ipsillateral & upper vis field extent
%             pss.(sn).stimCtr = [12.5, -8, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [30, 30];  %35 *[1 1];

%             pss.(sn).stimCtr = [18, -5, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [35,35];  %35 *[1 1];

%             pss.(sn).stimCtr = [18, -7.5, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [30,30];  %35 *[1 1];
            
%             pss.(sn).stimCtr = [11, -11, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [20,20];  %35 *[1 1];

            pss.(sn).stimCtr = [12, -12, 0]; % [x, y, zOffsetFromFixation]
            pss.(sn).gridSz  = [20, 20];  %35 *[1 1];

%             pss.(sn).stimCtr = [12, -10, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [20,20];  %35 *[1 1];
%             
%             pss.(sn).stimCtr = [14, -8, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [20, 40];  %35 *[1 1];
            
            
%             % fwhm==  3    @ 9x9
%             pss.(sn).stimCtr = [15, -8.3, 0]; % [x, y, zOffsetFromFixation]
%             pss.(sn).gridSz  = [35,30];  %35 *[1 1];
            
%             %             % HUGE searching grid
%                         pss.(sn).stimCtr = [20, 0, 0]; % [x, y, zOffsetFromFixation]
%                         pss.(sn).gridSz  = [35, 40];% *[1 1];
            
        end
        % zero out depth of grid for now...
        pss.(sn).stimCtr(3) = 0;%pss.display.viewdist;
        pss.(sn).gridSz(3) = 3;
        
        % Rendering Space:
        %   'enviro'    --3D translations & rotations remain in envirocentric coordinates
        %                 i.e. translation w/in display plane, and motion directions are relative to frontoparallel plane
        %   'direction' --Same spatial coords as 'enviro' (i.e. translation w/in display plane),
        %                 but motion directions are relative to subject 
        %   'retino'    --Translation & motion direction always relative to observer
        %                 "toward" always toward subject, and eccentric locations rendered with disparity
        %                 to produce stimulus viewing distance matched to fixation
        
        pss.(sn).renderingSpace = 'retino';
        
        % spherical dot volume
        pss.(sn).ballSz = max(pss.(sn).gridSz(1:2)./(gridN-1)) * 1.3;% % dot ball volume dia. in DEG   %min(pss.(sn).gridSz(1:2)./(gridN-1)) * sqrt(2);   %8;

        % Shared dot parameters ("snBase")
        pss.(sn).ndots = 120; % 200;%
        pss.(sn).dotType = 5;
        pss.(sn).dotSz = 0.4; % dot dia.:  0.3 if dotType<=2, or 0.5 CM if dotType>2) 
        pss.(sn).dotSpd = 8; % deg/sec

        pss.(sn).col = kron(randi([0,1], [1, pss.(sn).ndots]), [1 1 1]');
        pss.(sn).col(4,:) =  1; % all fully opaque

        % dot ball limits (types:  'cart')
        pss.(sn).limtype = 'cart';
                
        pss.(sn).lims = kron( pss.(sn).ballSz, 0.5*[-1 -1 -1; 1 1 1]); % [min; max] limits (cm)
        
        % overlay markers & mouse tracking
        pss.(sn).drawMarkers = true;
        pss.(sn).markerRects = [];
        pss.(sn).trackMouse = false; % no mouse tracking for RfPos stim
        
        % initial volume rotation so that motion in the "[0 0]" direction is rightward frontoparallel
        rz = d2r(90);
        initialRotMatrix = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];

        % *** Total stimulus duration [motionDur] set during pmBase module setup above ***
        
        tmpModule = pldapsModule('modName',sn, 'name','glDraw.dotBall_beta', 'matrixModule',true, 'order',10);
        
        % adjustments & defaults
        tmpModule.use = false;
        tmpModule.on = false;
        
        tmpModule.pos = zeros(3, pss.(sn).ndots);
        tmpModule.dir = zeros(2, 1);    %randi(180, 1, pss.(sn).ngabors);
        tmpModule.stimPos = [0, 0, 0]; % stim position relative to stimCtr
        tmpModule.frameRotMatrix = initialRotMatrix;

        % Timing
        tmpModule.isi = .0;
        
        ncopies = 12;
        motionModuleDur = fixDur/ncopies;
        
        % Report to command window
        fprintLineBreak;
        sr = 1000;
        fprintf('\t~~~\tmotionDur: %3.2fms,  isi: %3.2fms  ...x%d== %2.2fs total\n', (motionModuleDur-tmpModule.isi)*sr, tmpModule.isi*sr, ncopies, fixDur);
        fprintLineBreak;
        
        
    otherwise
        error('Unrecognized stimMode requested.')
end

% create duplicate/indexed stim modules for each presentation
matrixModNames = {};
for i = 1:ncopies
    mN = sprintf('%s%0.2d', sn, i);
    pss.(mN) = tmpModule;
    pss.(mN).stateFunction.modName = mN;
    pss.(mN).stateFunction.matrixModule = true;
    % timing: each module [onset, offset] time in sec (relative to MOTION state start)
    basedur = (i-1)*motionModuleDur;
    pss.(mN).modOnDur = [0, motionModuleDur-tmpModule.isi] +basedur; % subtract isi from offset
    % module names [for local use]
    matrixModNames{i} = mN;
end

%% Create Pldaps object structure
p = pldaps(subj, pss);


%% Generate conditions matrix

% Drift directions
dinc = 90; % radial location increment
% dirs = dinc:dinc:360; % radial locations
% dirs = [dirs, 90,270]; % append TOWARD & AWAY (when volume is tilted by 90 deg into xz axis)
% dirs = [0, 180, 90,270]; % append TOWARD & AWAY (when volume is tilted by 90 deg into xz axis)

dirs = (dinc:dinc:360)';
dirs = [0, 90]';
dirs(:,2) = 0;
% dirs = [dirs; [90,90; 270,90]]; % append TOWARD & AWAY (when volume is tilted by 90 deg into xz axis)

% Cartesian grid locations 
% -- center & span moved into matrixModule base
% grid0   = [-10, 10];    % grid center  (visual deg)
% gridSz  = 12 * [1 1];   % Full width of grid  (visual deg)
% number of samples in xy dimension
% gridN   = ceil(gridSz./pss.(matrixModNames{1}).gabFwhm) * [1 1];
xs = linspace(-.5,.5, gridN(1));     % .*gridSz(1)  +grid0(1);
ys = linspace(-.5,.5, gridN(2));     % .*gridSz(2)  +grid0(2);

xs(1:ipsiColumn(1)) = 1.3*xs(1:ipsiColumn(1));
ys(1:ipsiColumn(2)) = 1.3*ys(1:ipsiColumn(2));

if pss.(sn).gridSz(3) <= 1
    depthCM = 0
else
    depthCM = disp2depth( 60*linspace(-1,1, pss.(sn).gridSz(3))', pss.display.viewdist/100, pss.display.ipd)
end
% make a fully crossed matrix
[xx, yy, zz, dd] = ndgrid( xs, ys, depthCM, 1:size(dirs,1));

% is2d = numel(dd(:,:,1:end-2));% make last two dirs opposite binocular

%% Use new style condMatrix    
% name of pldaps module that contains experimental manipulation
% ...how could we better know/extract this when needed?
% sn = tstMode;

% Make the condMatrix!
c = cell(size(xx)); % *** maintain same shape as condition matrix source values

for i = 1:numel(c) % 1:length(dirs)
    % Set up conditions
    % stim position relative to center
    c{i}.stimPos    = [ xx(i), yy(i), zz(i)];
    % motion direction [direction, xz tilt]
    c{i}.dir        = dirs(dd(i),:);%* ones(tmpModule.ngabors, 2);

end

% Add condition matrix to the pldaps structure
%   p.conditions = c; % !! leave this empty, everything is in .condMatrix now
p.condMatrix.conditions = c;

%% Setup Block adjustment of viewdist
% vdists = [1, 2] .* p.trial.display.viewdist;

if exist('vdists','var') && ~isempty(vdists)
    % SAFETY CHECK
    % - TODO: build this into rig/code!
    if any(vdists<30) || any(vdists>120)
        error('ViewDist:Setup','Requested viewdist(s) %s exceeds dynamic range',mat2str(vdists));
    end
    
    B = cell(size(vdists));
    for i = 1:numel(vdists)
        B{i}.display.viewdist = vdists(i);
    end
    
    p.condMatrix.blocks = B;
end


%% Initialize condMatrix (w/ randMode, pass, & block parameters)
% p.condMatrix = condMatrix(p, 'randMode', [1,2,3,4], 'nPasses',inf, 'blockModulo',3);  % proper experiment mode
p.condMatrix = condMatrix(p, 'randMode', [1], 'nPasses',inf, 'blockModulo',2);    % DEBUG


%% Run it!!
p.run;



end
