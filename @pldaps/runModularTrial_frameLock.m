function p = runModularTrial_frameLock(p)
%runModularTrial    runs a single Trial by calling the functions defined in 
%            their substruct (found by pldaps.getModules) and the function
%            defined in p.trial.pldaps.trialFunction through different states
%
% 03/2013 jly     Wrote hyperflow
% 03/2014 jk      Used jly's code to get the PLDAPS structure and frame it into a class
%                 might change to ASYNC buffer flipping. but won't for now.
% 03/2016 jk      modular version
% 2019-xx-xx TBC  Separate version to artifically impose frame-synced timing so that
%                 slow operations (like movie capture) can occur without trial time elapsing
% 2020-11-05 TBC  Cleaned & commented (somewhat)
% 

    persistent beenWarned
    if isempty(beenWarned)
        fprintLineBreak('*!*',30)
        fprintLineBreak('*!*',30)
        fprintf(2, 'PLDAPS setup to use trial function:  %s\nThis will frame lock stimulus "time" to monitor refresh rate (for movie creation),\nand is not be suitable for experimental presentation.\nYou have been warned...\n', mfilename);
        fprintLineBreak('*!*',30)
        fprintLineBreak('*!*',30)
        beenWarned = 1;
    end
    
    %get all functionHandles that we want to use
    [modules, moduleFunctionHandles, moduleRequestedStates] = getModules(p);

    %order the framestates that we will iterate through each trial by their value
    % only positive states are frame states. And they will be called in
    % order of the value. Check comments in pldaps.getReorderedFrameStates
    % for explanations of the default states
    % negative states are special states outside of a frame (trial,
    % experiment, etc)
    [stateValue, stateName] = p.getReorderedFrameStates(p.trial.pldaps.trialStates, moduleRequestedStates);
    nStates=length(stateValue);

    %% trialSetup
    runStateforModules(p, 'trialSetup', modules, moduleFunctionHandles, moduleRequestedStates);

    %% trialPrepare
    %   called just before the trial starts for time critical calls
    %   (e.g. to start data aquisition)
    runStateforModules(p, 'trialPrepare', modules, moduleFunctionHandles, moduleRequestedStates);

    %%% MAIN WHILE LOOP %%%
    %-------------------------------------------------------------------------%
    while ~p.trial.flagNextTrial && p.trial.pldaps.quit == 0
        % Advance to next frame, update frame index
        p.trial.iFrame = p.trial.iFrame + 1;
        
        % time of the next flip request
        %   --see pldapsDefaultTrial.m>>frameFlip to determine how/if this gets applied
        p.trial.nextFrameTime = p.trial.timing.timeLastFrame + p.trial.display.ifi;

        % Start timer for GPU rendering operations
        Screen('GetWindowInfo', p.trial.display.ptr, 5);
        
        %% PLDAPS frame states
        % iterate through [.frameUpdate, .framePrepareDrawing, .frameDraw, ...etc]
        for iState=1:nStates
            runStateforModules(p, stateName{iState}, modules, moduleFunctionHandles, moduleRequestedStates);

            %% FUDGE the trial time as if only one frame has elapsed regardless of how much
            % processing/pause time has elapsed
            % (important for recording stimulus movie frames or debugging code)
            p.trial.ttime = p.trial.iFrame*p.trial.display.ifi + p.trial.trstart;
            
            p.trial.remainingFrameTime = p.trial.nextFrameTime - p.trial.ttime;
            p.trial.timing.frameStateChangeTimes(iState, p.trial.iFrame) = p.trial.ttime - p.trial.nextFrameTime + p.trial.display.ifi;
        end
        
        % Retrieve GPU render time of last frame
        dinfo = Screen('GetWindowInfo', p.trial.display.ptr);
        p.trial.frameRenderTime(p.trial.iFrame) = dinfo.GPULastFrameRenderTime;

    end %while Trial running

    %% trialItiDraw
    %  ** Inherently not a time-critical operation, so no call to setTimeAndFrameState necessary
    %   ...also, setTimeAndFrameState uses current state as an index, so using with this would break
    runStateforModules(p, 'trialItiDraw', modules, moduleFunctionHandles, moduleRequestedStates);

    if round(Priority)<MaxPriority('GetSecs')
        warning('pldaps:runTrial', 'Thread priority was degraded by operating system during the trial.')
    end

    %% trialCleanUpandSave
    runStateforModules(p, 'trialCleanUpandSave', modules, moduleFunctionHandles, moduleRequestedStates);

end %runModularTrial
