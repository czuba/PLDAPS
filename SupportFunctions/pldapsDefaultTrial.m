function pldapsDefaultTrial(p, state, sn)

    if nargin<3
        error('Use of PLDAPS modules without 3rd input string for module name is no longer allowed. See also  pldapsModule.m')
        % default p.trial.stimulus field was killed off a while back --TBC 12/2019
        %         sn='stimulus';
    end
    
    switch state
        % FRAME STATES
        % - each of these states are invoked on every display frame flip
        % (always place these at the top of your switch statement b/c they happen most frequently)
        case p.trial.pldaps.trialStates.frameUpdate
            frameUpdate(p);
            
            
        case p.trial.pldaps.trialStates.framePrepareDrawing 
            % ...does nothing by default
            %    framePrepareDrawing(p);
            %
            % This is where you should compute your stimulus features for the next frame.
            % At this point, the eye position, datapixx info, spike signals (maybe more)
            % should be as up-to-date as they are going to be before the next display frame
            %
            % Currently empty in the pldapsDefaultTrialFunction by design:
            % This is where YOUR code should do most of the operations necessary for
            % each display frame. E.g.
            % - assess eye position relative to [fixation] elements onscreen
            % - set/update any parameters for the upcoming stimulus frame
            % - determine whether behavioral conditions are met,
            %   and/or the trial has been completed
            %
            
            
        case p.trial.pldaps.trialStates.frameDraw
            for i = p.trial.display.bufferIdx
                frameDraw(p, sn, i);
            end
            
            
        case p.trial.pldaps.trialStates.frameGLDrawLeft
            % Skip altogether if not enabled
            if p.trial.display.useGL
                % Only pass in the needed display structure
                frameGLDrawLeft(p.trial.display);
            end
            
            
        case p.trial.pldaps.trialStates.frameGLDrawRight
            % Skip altogether if not enabled & not stereomode
            if p.trial.display.useGL && p.trial.display.stereoMode~=0
                % Only pass in the needed display structure
                frameGLDrawRight(p.trial.display);
            end
            
            
        case p.trial.pldaps.trialStates.frameDrawingFinished
            frameDrawingFinished(p);
            
            
        case p.trial.pldaps.trialStates.frameFlip
            frameFlip(p);
            
            
        % TRIAL STATES
        % - each invoked once per trial, at the corresponding start/end of every experimental trial
        case p.trial.pldaps.trialStates.trialItiDraw
            trialItiDraw(p);
            % This state allows drawing to the screen during the intertrial interval(iti),
            % just prior to the final screen flip that occurs during trialCleanUpandSave.
            % Anything drawn to the screen during this state will remain [static] on screen between trials.
            % - good for things like fiducials or static UI elements that can be distracting
            %   if they blink out of existance between trials
            % ** Note: use this state with caution, so that the subject/animal isn't left
            %    unintentionally fixating a target without any ongoing tracking or reward potential!
            
            
        case p.trial.pldaps.trialStates.trialSetup
            trialSetup(p);
            
            
        case p.trial.pldaps.trialStates.trialPrepare
            trialPrepare(p);
            
            
        case p.trial.pldaps.trialStates.trialCleanUpandSave
            cleanUpandSave(p);
            
            
        % EXPERIMENT STATES
        % - only invoked once at start/end of PLDAPS experiment session
        case p.trial.pldaps.trialStates.experimentPreOpenScreen
            experimentPreOpenScreen(p, sn)
            
            
        case p.trial.pldaps.trialStates.experimentAfterTrials
            if ~isempty(p.trial.pldaps.experimentAfterTrialsFunction)
               h=str2func(p.trial.pldaps.experimentAfterTrialsFunction);
               h(p, state);
            end
            
            
        case p.trial.pldaps.trialStates.experimentPostOpenScreen
            experimentPostOpenScreen(p, sn);
            
                
    end
    
end %main function


% % % % % % % % % % % % % % % 
% % % Sub-functions
% % % % % % % % % % % % % % % 

%---------------------------------------------------------------------% 
%% frameUpdate   (check & refresh keyboard/mouse/analog/eye data)
function frameUpdate(p)
%%TODO: add buffer for Keyboard presses, mouse position and clicks.

% Check keyboard
p = pds.keyboard.getQueue(p);

% Some standard PLDAPS key functions
if any(p.trial.keyboard.firstPressQ)
    
    % [M]anual reward
    if  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.mKey)
        pds.behavior.reward.give(p);
        
        % [P]ause
    elseif  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.pKey)
        p.trial.pldaps.quit = 1;
        ShowCursor;
        
        % [Q]uit
    elseif  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.qKey)
        p.trial.pldaps.quit = 2;
        ShowCursor;
        
        % [D]ebug mode   (...like pause, but does not leave workspace of currently executing trial)
    elseif  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.dKey)
        disp('stepped into debugger. Type return to start first trial...')
        keyboard %#ok<MCKBD>
        
        % [C]alibration
        % NOTE:  This OTF calibration trigger did not work in practice
        % -- Necessary to start a new trial for calibration & changes to parameters here
        % don't readily carry over to subsequent trials(by PLDAPS design)
        % -- Should also require a modifier key to reduce chance of inadvertently triggering
        %
        % % %             elseif  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.cKey)
        % % %                 % data interruption/corruption
        % % %                 p.trial.flagNextTrial = true;
        % % %                 p.trial.tracking.on = true;
        
    end
end

% Poll mouse
if p.trial.mouse.use
    [cursorX,cursorY,buttonState] = GetMouse(p.trial.mouse.windowPtr);
    % Return data in trial struct
    p.trial.mouse.samples = p.trial.mouse.samples+1;
    p.trial.mouse.samplesTimes(p.trial.mouse.samples)=GetSecs;
    p.trial.mouse.cursorSamples(1:2,p.trial.mouse.samples) = [cursorX;cursorY];
    p.trial.mouse.buttonPressSamples(:,p.trial.mouse.samples) = buttonState';
    
    % Use as eyepos if requested
    if p.trial.mouse.useAsEyepos
        if p.trial.pldaps.eyeposMovAv==1
            p.trial.eyeX = p.trial.mouse.cursorSamples(1,p.trial.mouse.samples);
            p.trial.eyeY = p.trial.mouse.cursorSamples(2,p.trial.mouse.samples);
        else
            mInds=(p.trial.mouse.samples-p.trial.pldaps.eyeposMovAv+1):p.trial.mouse.samples;
            p.trial.eyeX = mean(p.trial.mouse.cursorSamples(1,mInds));
            p.trial.eyeY = mean(p.trial.mouse.cursorSamples(2,mInds));
        end
        % Also report delta "eye" position
        nback = p.trial.mouse.samples + [-1,0]; nback(nback<1) = 1;
        p.trial.eyeDelta = [diff(p.trial.mouse.cursorSamples(1,nback), [], 2),...
            diff(p.trial.mouse.cursorSamples(2,nback), [], 2)];
        
    end
end

%get analogData from Datapixx
pds.datapixx.adc.getData(p);

%get eyelink data
if p.trial.eyelink.use      && ~p.trial.tracking.use
    pds.eyelink.getQueue(p);
end

if p.trial.tracking.use
    % update from source & apply calibration
    pds.tracking.frameUpdate(p);
end
%get plexon spikes
% pds.plexon.spikeserver.getSpikes(p);

end %frameUpdate


%---------------------------------------------------------------------%
%%  framePrepareDrawing
function framePrepareDrawing(p) %#ok<INUSD,DEFNU>
% This is where you should compute your stimulus features for the next frame.
% At this point, the eye position, datapixx info, spike signals (maybe more)
% should be as up-to-date as they are going to be before the next display frame
%
% Currently empty in the pldapsDefaultTrialFunction by design:
% This is where YOUR code should do most of the operations necessary for
% each display frame. E.g.
% - assess eye position relative to [fixation] elements onscreen
% - set/update any parameters for the upcoming stimulus frame
% - determine whether behavioral conditions are met,
%   and/or the trial has been completed
%

end %framePrepareDrawing


%---------------------------------------------------------------------%
%%  frameDraw
function frameDraw(p, sn, drawBuffer)
%this holds the code to draw some stuff to the overlay (using
%switches, like the grid, the eye Position, etc

% be smart(ish) about binocular rendering
Screen('SelectStereoDrawBuffer', p.trial.display.ptr, drawBuffer);

% Grid overlay
if p.trial.pldaps.draw.grid.use
    Screen('DrawLines', p.trial.display.overlayptr, p.trial.pldaps.draw.grid.tick_line_matrix, 1, p.trial.display.clut.window, p.trial.display.ctr(1:2));
end

% Framerate history (only render once, to left eye)
if p.trial.pldaps.draw.framerate.use && p.trial.iFrame>2 && ~drawBuffer
    pds.pldaps.draw.framerateDraw(p)
end

% Eye positon
if  p.trial.pldaps.draw.eyepos.use
    if ~p.trial.eyelink.use || numel(p.trial.eyelink.eyeIdx)==1
        Screen('Drawdots', p.trial.display.overlayptr, [p.trial.eyeX p.trial.eyeY]', ...
            p.trial.(sn).eyeW, p.trial.display.clut.eyepos, [0 0],0);
    else
        Screen('Drawdots', p.trial.display.overlayptr, [p.trial.eyeX(drawBuffer+1) p.trial.eyeY(drawBuffer+1)]', ...
            p.trial.(sn).eyeW, p.trial.display.clut.(['eye',num2str(drawBuffer)]), [0 0],0);
    end
end

if p.trial.mouse.use && p.trial.pldaps.draw.cursor.use && ~drawBuffer
    Screen('Drawdots',  p.trial.display.overlayptr,  p.trial.mouse.cursorSamples(1:2,p.trial.mouse.samples), ...
        p.trial.(sn).eyeW, p.trial.display.clut.cursor, [0 0],0);
end

% Photodiode sync flash
if p.trial.pldaps.draw.photodiode.use && mod(p.trial.iFrame, p.trial.pldaps.draw.photodiode.everyXFrames) == 0
    p.trial.timing.photodiodeTimes(:,p.trial.pldaps.draw.photodiode.dataEnd) = [p.trial.ttime p.trial.iFrame];
    p.trial.pldaps.draw.photodiode.dataEnd = p.trial.pldaps.draw.photodiode.dataEnd+1;
    Screen('FillRect', p.trial.display.ptr, [1 1 1]', p.trial.pldaps.draw.photodiode.rect');
end
end %frameDraw


%---------------------------------------------------------------------%
%%  frameGLDrawLeft
function frameGLDrawLeft(ds)
global GL

% Switch to 3D mode:
Screen('BeginOpenGL', ds.ptr);

% Select stereo draw buffer WITHIN a 3D openGL context!
% unbind current FBOS first     (per PTB source:  "otherwise bad things can happen...")
glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, uint32(0))

fbo = uint32(1); % FBO index for RIGHT eye
glBindFramebufferEXT(GL.READ_FRAMEBUFFER_EXT, fbo);
glBindFramebufferEXT(GL.DRAW_FRAMEBUFFER_EXT, fbo);
glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, fbo);

% Clear depth (and color?) buffers:
glClear(GL.DEPTH_BUFFER_BIT);

% Setup camera for this eyes 'view':
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
% OBSERVER-CENTRIC; (0,0,0) should be observer location
%   gluLookAt( eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ )
%   !!**!! OpenGL space z-space is in the negative direction (else z-buffer doesn't work!)
gluLookAt(ds.obsPos(1)-0.5*ds.ipd*cosd(ds.obsPos(4)), ds.obsPos(2), -(ds.obsPos(3)-.5*ds.ipd*sind(ds.obsPos(4))),...
    ds.fixPos(1), ds.fixPos(2), -ds.fixPos(3),...
    ds.upVect(1), ds.upVect(2), ds.upVect(3));

end%  frameGLDrawLeft


%---------------------------------------------------------------------%
%%  frameGLDrawRight
function frameGLDrawRight(ds)
global GL

% We are already inside 3D openGL 'BeginOpenGL' mode at this point
% Select stereo draw buffer WITHIN a 3D openGL context!
% unbind current FBOS first     (per PTB source:  "otherwise bad things can happen...")
glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, uint32(0))

fbo = uint32(2); % FBO index for RIGHT eye
glBindFramebufferEXT(GL.READ_FRAMEBUFFER_EXT, fbo);
glBindFramebufferEXT(GL.DRAW_FRAMEBUFFER_EXT, fbo);
glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, fbo);

% Clear depth buffers:  (CONFIRMED:  This must be done for each eye on each frame. --TBC)
glClear(GL.DEPTH_BUFFER_BIT);

% Setup camera for RIGHT eye 'view':
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
% OBSERVER-CENTRIC; (0,0,0) should be observer location
%   gluLookAt( eyeX, eyeY, -eyeZ, centerX, centerY, -centerZ, upX, upY, upZ )
%   !!**!! OpenGL space z-space is in the negative direction (else z-buffer doesn't work!)
gluLookAt(ds.obsPos(1)+0.5*ds.ipd*cosd(ds.obsPos(4)), ds.obsPos(2), -(ds.obsPos(3)+.5*ds.ipd*sind(ds.obsPos(4))),...
    ds.fixPos(1), ds.fixPos(2), -ds.fixPos(3),...
    ds.upVect(1), ds.upVect(2), ds.upVect(3));

end %frameGLDrawRight


%---------------------------------------------------------------------%
%%  frameDrawingFinished
function frameDrawingFinished(p)
% Check for & disable OpenGL mode before any Screen() calls (else will crash)
if p.trial.display.useGL
    % [~, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');
    % % NOTE: polling opengl for info is slow, so make sense to go-for-broke here
    Screen('EndOpenGL', p.trial.display.ptr);
end

Screen('DrawingFinished', p.trial.display.ptr);

if p.trial.datapixx.use && ~isempty(p.trial.datapixx.strobeQ)
    % send all the pending strobes and clear the queue
    pds.datapixx.strobeQueue(p.trial.datapixx.strobeQ);
    p.trial.datapixx.strobeQ = [];
end

end %frameDrawingFinished


%---------------------------------------------------------------------%
%%  frameFlip
function frameFlip(p)
ft=cell(5,1); % empty cell to receive all outputs from Screen('Flip'...)
% nextFrameTime = p.trial.nextFrameTime+p.trial.trstart -0.05*p.trial.display.ifi;
nextFrameTime = 0; % flip on next available frame

[ft{:}] = Screen('Flip', p.trial.display.ptr, nextFrameTime);
p.trial.timing.flipTimes(:,p.trial.iFrame) = [ft{:}];

% The overlay screen always needs to be initialized with a FillRect call
if p.trial.display.overlayptr ~= p.trial.display.ptr
    Screen('FillRect', p.trial.display.overlayptr, 0);
end

p.trial.timing.timeLastFrame = p.trial.timing.flipTimes(1,p.trial.iFrame)-p.trial.trstart;

end %frameFlip


%---------------------------------------------------------------------%
%%  trialSetup
function trialSetup(p)
% ensure max priority at onset of every trial
% Switch to high priority mode
if p.trial.pldaps.maxPriority
    % just set it already!
    if IsLinux
        Priority(p.trial.pldaps.maxPriority);
    else
        Priority(MaxPriority('GetSecs'));
    end

    %     oldPriority=Priority;
    %     maxPriority=MaxPriority('GetSecs');
    %     if oldPriority < maxPriority
    %         Priority(maxPriority);
    %     end
end

% fprintf('\n');  %DEBUG

% p.trial.display = p.static.display.syncToTrialStruct(p.trial.display);
% disp(p.trial.display.fixPos)
% disp(p.trial.display.obsPos)
p.trial.timing.flipTimes       = zeros(5,p.trial.pldaps.maxFrames);
p.trial.timing.frameStateChangeTimes=nan(9,p.trial.pldaps.maxFrames);

if(p.trial.pldaps.draw.photodiode.use)
    p.trial.timing.photodiodeTimes=nan(2,p.trial.pldaps.maxFrames);
    p.trial.pldaps.draw.photodiode.dataEnd=1;
end

%these are things that are specific to subunits as eyelink,
%datapixx, mouse and should probabbly be in separarte functions,
%but I have no logic/structure for that atm.

if p.trial.datapixx.use
    %setup analogData collection from Datapixx
    pds.datapixx.adc.trialSetup(p);
    % Sync Datapixx & PTB clocks (via streamlined version of PsychDataPixx('GetPreciseTime'))
    p.trial.timing.datapixxPreciseTime = pds.datapixx.syncClocks(p.trial.datapixx.GetPreciseTime);
end

% setup  fields for keyboard data
pds.keyboard.trialSetup(p);

% setup tracking calibration
if p.trial.tracking.use
    pds.tracking.trialSetup(p);
end

%setup a fields for the mouse data
if p.trial.mouse.use
    pds.mouse.trialSetup(p);
end

%%% Eyelink Toolbox Setup %%%
%-------------------------------------------------------------------------%
% preallocate for all eye samples and event data from the eyelink
if p.trial.eyelink.use && ~p.trial.tracking.use
    pds.eyelink.startTrial(p);
end

%%% Spike server
%-------------------------------------------------------------------------%
if  p.trial.plexon.spikeserver.use
    % NOTE: This functionality/setup has been lost in the sands of time (circa Oct. 2020)
    p.trial.plexon.spikeserver.spikeCount=0;
    pds.plexon.spikeserver.getSpikes(p); %save all spikes that arrives in the inter trial interval
end

%%% prepare reward system
pds.behavior.reward.trialSetup(p);

%%% prepare to plot framerate history onscreen
if p.trial.pldaps.draw.framerate.use
    % Legacy debug behavior [YMMV]
    pds.pldaps.draw.framerateTrialSetup(p);
end

if p.trial.display.useGL
    % MOVED this step to common code:   pds.display.setViewdist.m
    %       setupGLPerspective(p.trial.display); % subfunction
    pds.display.setViewdist(p);
end

end %trialSetup


%---------------------------------------------------------------------%
%%  trialPrepare
function trialPrepare(p)

%%% setup PsychPortAudio %%%
%-------------------------------------------------------------------------%
pds.sound.clearBuffer(p)

% Ensure anything in the datapixx buffer has been pushed/updated
if p.trial.datapixx.use
    Datapixx RegWrRd;
end

%%% Initalize Keyboard %%%
%-------------------------------------------------------------------------%
pds.keyboard.clearBuffer(p);

%%% Eyelink Toolbox Setup %%%
%-------------------------------------------------------------------------%
% preallocate for all eye samples and event data from the eyelink
pds.eyelink.startTrialPrepare(p);


%%% START OF TRIAL TIMING %%
%-------------------------------------------------------------------------%
% record start of trial in Datapixx, PLDAPS & Plexon
% each device has a separate clock

% At the beginning of each trial, strobe a unique number to plexon
% through the Datapixx to identify each trial. Often the Stimulus display
% will be running for many trials before the recording begins so this lets
% the plexon rig sync up its first trial with whatever trial number is on
% for stimulus display.
% SYNC clocks

% Construct a unique trial number based on:
%   6 element clock time: [year, month, day, 24hour, minute, second]
unique_number = fix(clock);
% 	substitute year with trial number (i.e. something actually relevant on the scale of an experimental session)
unique_number(1) = p.trial.pldaps.iTrial;
%   shift unique numbers into the upper half of our 15-bit strobed word range
unique_number = unique_number + 2^14;
%  ...leaving the lower 16,383 values for easily identifiable event values (p.trial.event).

if p.trial.datapixx.use
    % Strobe unique_number via datapixx
    for i = 1:numel(unique_number)
        p.trial.datapixx.unique_number_time(i,:) = pds.datapixx.strobe(unique_number(i));
    end
end
p.trial.unique_number = unique_number;    % trial identifier

if p.trial.datapixx.use
    % start of trial sync signal (Plexon)
    p.trial.timing.datapixxStartTime = Datapixx('Gettime');
    p.trial.timing.datapixxTRIALSTART = pds.datapixx.strobe(p.trial.event.TRIALSTART);
end


% These params are all predetermined, so just set them equal to 0,
% and keep any code post-vblsync to an absolute minimum!  (...yes, even just touching p.trial)
p.trial.ttime  = 0;                     % formerly:  GetSecs - p.trial.trstart;
p.trial.timing.syncTimeDuration = 0;    % formerly:  p.trial.ttime;
p.trial.timing.timeLastFrame = 0;       % formerly:  vblTime-p.trial.trstart;

% Sync up with screen refresh before jumping into actual trial
%   ** this also ensures that the async flip scheduled at the end of the last trial
%      has had time to complete & won't interfere with future draws/flips
p.trial.timing.itiFrameCount = Screen('WaitBlanking', p.trial.display.ptr);
p.trial.trstart = GetSecs;

% Tell datapixx to save a timestamp marker at the next frame flip. These will be
% transferred from datapixx box to PTB machine during completion of experiment (by run.m)
PsychDataPixx('LogOnsetTimestamps',1);%2

end %trialPrepare


%---------------------------------------------------------------------%
%%  trialItiDraw
function p = trialItiDraw(p)
% Only do the basic drawing commands here
%   ...e.g. maybe not eye pos, since it cannot be updated during this phase

% Be smart(ish) about binocular rendering
Screen('SelectStereoDrawBuffer', p.trial.display.ptr, 1);

% Grid overlay
if p.trial.pldaps.draw.grid.use
    % render grid to monocular or left eye buffer
    Screen('DrawLines', p.trial.display.overlayptr, p.trial.pldaps.draw.grid.tick_line_matrix, 1, p.trial.display.clut.window, p.trial.display.ctr(1:2));
    if any( (2:5)==p.trial.display.stereoMode )
        % also render overlay grid to other eye for non-overlapping stereomodes
        Screen('SelectStereoDrawBuffer', p.trial.display.ptr, 1);
        Screen('DrawLines', p.trial.display.overlayptr, p.trial.pldaps.draw.grid.tick_line_matrix, 1, p.trial.display.clut.window, p.trial.display.ctr(1:2));
    end
end

end %trialItiDraw


%---------------------------------------------------------------------%
%%  cleanUpandSave
function p = cleanUpandSave(p)

% Schedule a flip to occur at the next possible time, but don't bother waiting around for it.
%         Screen('AsyncFlipBegin', p.trial.display.ptr);
Screen('Flip', p.trial.display.ptr, 0, [], 1);
% Whatever was drawn to this screen will be visible throughout the inter-trial interval.
% This was previously always/only a blank screen, but if you want anything present on
% this screen, it can now be drawn during the  .trialItiDraw  state.
% NOTE: This is not a time-critical draw, and async flips do not return a valid timestamp
%       at time of schedule.

% Wait a fraction of a sec for any pending strobes to register downstream
WaitSecs(1e-4);     % 'UntilTime', t0+3e-5);

% Execute all time-sesitive tasks first
if p.trial.datapixx.use
    p.trial.datapixx.datapixxstoptime = Datapixx('GetTime');
end
p.trial.trialend = GetSecs- p.trial.trstart;

%clean up analogData collection from Datapixx
pds.datapixx.adc.cleanUpandSave(p);
if p.trial.datapixx.use
    % end of trial sync signal (Plexon)
    p.trial.timing.datapixxTRIALEND = pds.datapixx.strobe(p.trial.event.TRIALEND);
end

if(p.trial.pldaps.draw.photodiode.use)
    p.trial.timing.photodiodeTimes(:,p.trial.pldaps.draw.photodiode.dataEnd:end)=[];
end


p.trial.trialnumber   = p.trial.pldaps.iTrial;
p.trial.timing.flipTimes      = p.trial.timing.flipTimes(:,1:p.trial.iFrame);
p.trial.timing.frameStateChangeTimes    = p.trial.timing.frameStateChangeTimes(:,1:p.trial.iFrame);

%do a last frameUpdate   (checks & refreshes keyboard/mouse/analog/eye data)
frameUpdate(p)

%         % Flush KbQueue
%         KbQueueStop();
%         KbQueueFlush();

%will this crash when more samples where created than preallocated?
% mouse input
if p.trial.mouse.use
    i0 = p.trial.mouse.samples+1;
    p.trial.mouse.cursorSamples(:,i0:end) = [];
    p.trial.mouse.buttonPressSamples(:,i0:end) = [];
    p.trial.mouse.samplesTimes(:,i0:end) = [];
end

i0 = p.trial.keyboard.samples+1;
p.trial.keyboard.samplesTimes(:,i0:end) = [];
p.trial.keyboard.samplesFrames(:,i0:end) = [];
p.trial.keyboard.pressedSamples(:,i0:end) = [];
p.trial.keyboard.firstPressSamples(:,i0:end) = [];
p.trial.keyboard.firstReleaseSamples(:,i0:end) = [];
p.trial.keyboard.lastPressSamples(:,i0:end) = [];
p.trial.keyboard.lastReleaseSamples(:,i0:end) = [];


%---------------------------------------------------------------------%
% Plexon specific:
% Get spike server spikes
if p.trial.plexon.spikeserver.use
    try
        pds.plexon.spikeserver.getSpikes(p);
    catch me
        disp(me.message)
    end
end


%---------------------------------------------------------------------%
% Eyelink specific:
if p.trial.eyelink.use      && ~p.trial.tracking.use
    [Q, rowId] = pds.eyelink.saveQueue(p);
    p.trial.eyelink.samples = Q;
    p.trial.eyelink.sampleIds = rowId; % I overwrite everytime because PDStrialTemps get saved after every trial if we for some unforseen reason ever need this for each trial
    p.trial.eyelink.events   = p.trial.eyelink.events(:,~isnan(p.trial.eyelink.events(1,:)));
end


%---------------------------------------------------------------------%
% Reward system
pds.behavior.reward.cleanUpandSave(p);


end %cleanUpandSave


%% experimentPreOpenScreen(p)
function experimentPreOpenScreen(p, sn)
% setup some default modules

% ensure Priority setting is executable value
if p.trial.pldaps.maxPriority
    p.trial.pldaps.maxPriority = double(p.trial.pldaps.maxPriority);
    if IsLinux
        % maximum warp!!!
        p.trial.pldaps.maxPriority = 19;
    end
end
% Tracking module
if p.trial.tracking.use
    pds.tracking.setup(p, sn);
end

end %experimentPreOpenScreen


%% experimentPostOpenScreen(p)
function experimentPostOpenScreen(p, sn)
% Update screen-dependent params
if ~isfield(p.trial.(sn), 'eyeW')
    p.trial.(sn).eyeW = 8;
end

if ~isField(p.trial, 'event')
    defaultBitNames(p);
end

if p.trial.display.useGL
    % Initialize OpenGL rendering settings, viewport, & glPerspective
    % - force update for initialization
    pds.display.setViewdist(p, p.static.display.viewdist, 1);
end

% update p.trial with pdsDisplay object
p.trial.display = p.static.display.syncToTrialStruct(p.trial.display);

% % These are now dependent variables within pdsDisplay object (p.static.display)
% %         if isfield(p.trial.display, 'useGL') && p.trial.display.useGL
% %             %   .display.glPerspective == [fovy, aspect, zNear, zFar]
% %             p.trial.display.glPerspective = [atand(p.trial.display.heightcm/2/p.trial.display.viewdist)*2,...
% %                 p.trial.display.widthcm/p.trial.display.heightcm,...
% %                 p.trial.display.zNear,... % near clipping plane (cm)
% %                 p.trial.display.zFar];  % far clipping plane (cm)
% %         end


%---------------------------------------------------------------------%
% Identify tracking source & setup calibration fields
if p.trial.tracking.use
    pds.tracking.postOpenScreen(p);
end


end %experimentPostOpenScreen


% % %     %% setupGLPerspective
% % %     function setupGLPerspective(ds)
% % %         global GL
% % %
% % %         % readibility & avoid digging into this struct over & over
% % %         glP = ds.glPerspective;
% % %
% % %         % Setup projection matrix for each eye
% % %         % (** this does not change per-eye & [unlikely] between frames, so just do it once here)
% % %         Screen('BeginOpenGL', ds.ptr)
% % %
% % %         for view = 0%:double(ds.stereoMode>0)
% % %             % All of these settings will apply to BOTH eyes once implimented
% % %
% % %             % Select stereo draw buffer WITHIN a 3D openGL context!
% % %             % unbind current FBOS first (per PTB source:  "otherwise bad things can happen...")
% % %             glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, uint32(0))
% % %
% % %             % Bind this view's buffers
% % %             fbo = uint32(view+1);
% % %             glBindFramebufferEXT(GL.READ_FRAMEBUFFER_EXT, fbo);
% % %             glBindFramebufferEXT(GL.DRAW_FRAMEBUFFER_EXT, fbo);
% % %             glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, fbo);
% % %
% % %
% % %             % Setup projection for stereo viewing
% % %             glMatrixMode(GL.PROJECTION)
% % %             glLoadIdentity;
% % %             % glPerspective inputs: ( fovy, aspect, zNear, zFar )
% % %             gluPerspective(glP(1), glP(2), glP(3), glP(4));
% % %
% % %             % Enable proper occlusion handling via depth tests:
% % %             glEnable(GL.DEPTH_TEST);
% % %
% % %             % Enable alpha-blending for smooth dot drawing:
% % %             glEnable(GL.BLEND);
% % %             glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);
% % %
% % %             % 3D anti-aliasing?
% % %             % NOTE: None of the standard smoothing enables improve rendering of gluDisk or sphere.
% % %             % Only opening PTB window with multisampling (==4) has any effect
% % %
% % %             % basic colors
% % %             glClearColor(ds.bgColor(1), ds.bgColor(2), ds.bgColor(3), 1);
% % %             glColor4f(1,1,1,1);
% % %
% % %             % Disable lighting
% % %             glDisable(GL.LIGHTING);
% % %             % glDisable(GL.BLEND);
% % %
% % %             % % %             if ds.goNuts
% % %             % % %                 % ...or DO ALL THE THINGS!!!!
% % %             % % %                 % Enable lighting
% % %             % % %                 glEnable(GL.LIGHTING);
% % %             % % %                 glEnable(GL.LIGHT0);
% % %             % % %                 % Set light position:
% % %             % % %                 glLightfv(GL.LIGHT0,GL.POSITION, [1 2 3 0]);
% % %             % % %                 % Enable material colors based on glColorfv()
% % %             % % %                 glEnable(GL.COLOR_MATERIAL);
% % %             % % %                 glColorMaterial(GL.FRONT_AND_BACK, GL.AMBIENT_AND_DIFFUSE);
% % %             % % %                 glMaterialf(GL.FRONT_AND_BACK, GL.SHININESS, 48);
% % %             % % %                 glMaterialfv(GL.FRONT_AND_BACK, GL.SPECULAR, [.8 .8 .8 1]);
% % %             % % %             end
% % %         end
% % %         Screen('EndOpenGL', ds.ptr)
% % %     end %setupGLPerspective

