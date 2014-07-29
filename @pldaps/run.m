function dv = run(dv)
% [dv] = run()
% PLDAPS (Plexon Datapixx PsychToolbox) version 4
%       run is a wrapper for calling PLDAPS condition files
%           It opens PsychImaging pipeline and initializes datapixx for
%           dual color lookup tables. Everything else must be in the
%           condition file and trial function. See PLDAPScheck.m for
%           explanation.  IMPORTANT: edit setupPLDAPSenv.m and
%           makeRigConfigFile.m before running. 
% INPUTS:
%       subj [string]       - initials for subject
%       condition [string]  - name of matlab function that runs trials
%                           - you must have the condition file in your path 
%       newsession [0 or 1] - if 1, start new PDS. 0 load old PDS (defaults
%       to 1)
%       [brackets] indicate optional variables and their default values

% 10/2011 jly wrote it (modified from letsgorun.m)
% 12/2013 jly reboot. updated to version 3 format.
% 04/2014 jk  movied into a pldaps classed; adapated to new class structure

% Tested to run with Psychtoolbox
% 3.0.11 - Flavor: beta - Corresponds to SVN Revision 4331
% For more info visit:
% https://github.com/Psychtoolbox-3/Psychtoolbox-3

%TODO: 
% shoudl the outputfile uigetfile be optional?
% one same system for modules, e.g. moduleSetup, moduleUpdate, moduleClose
% make HideCursor optional
% TODO:reset class at end of experiment or mark as recorded, so I don't
% run the same again by mistake
% Todo save: defaultparameters beofre 1st trial
%done
% -change edit setupPLDAPSenv.m and  makeRigConfigFile.m before running.
% it's now createRigPrefs
% -make wait for return optional?: pldaps.pause.preExperiment=false

try
    %% Setup and File management
    % Enure we have an experimentSetupFile set and verify output file
    
    
    % pick YOUR experiment's main CONDITION file-- this is where all
    % expt-specific stuff emerges from
    if isempty(dv.defaultParameters.session.experimentSetupFile)
        [cfile, cpath] = uigetfile('*.m', 'choose condition file', [base '/CONDITION/debugcondition.m']); %#ok<NASGU>
        
        dotm = strfind(cfile, '.m');
        if ~isempty(dotm)
            cfile(dotm:end) = [];
        end
        dv.defaultParameters.session.experimentSetupFile = cfile;
    end
             
    dv.defaultParameters.session.initTime=now;
        
        
    if ~dv.defaultParameters.pldaps.nosave
        dv.defaultParameters.session.file = fullfile(dv.defaultParameters.pldaps.dirs.data, [dv.defaultParameters.session.subject datestr(dv.defaultParameters.session.initTime, 'yyyymmdd') dv.defaultParameters.session.experimentSetupFile datestr(dv.defaultParameters.session.initTime, 'HHMM') '.PDS']);
        [cfile, cdir] = uiputfile('.PDS', 'initialize experiment file', dv.defaultParameters.session.file);
        dv.defaultParameters.session.dir = cdir;
        dv.defaultParameters.session.file = cfile;
    else
        dv.defaultParameters.session.file='';
        dv.defaultParameters.session.dir='';
    end
        
    %% Open PLDAPS windows
    % Open PsychToolbox Screen
    dv = openScreen(dv);
    
    % Setup PLDAPS experiment condition
    dv = feval(dv.defaultParameters.session.experimentSetupFile, dv);
    
            %
            % Setup Photodiode stimuli
            %-------------------------------------------------------------------------%
            if(dv.trial.pldaps.draw.photodiode.use)
                makePhotodiodeRect(dv);
            end
    
            % Tick Marks
            %-------------------------------------------------------------------------%
            if(dv.trial.pldaps.draw.grid.use)
                dv = initTicks(dv);
            end


            %get and store changes of current code to the git repository
            dv = pds.git.setup(dv);
            
            %things that were in the conditionFile
            dv = pds.eyelink.setup(dv);
    
            %things that where in the default Trial Structure
            
            % Audio
            %-------------------------------------------------------------------------%
            dv = pds.audio.setup(dv);
            
            % Audio
            %-------------------------------------------------------------------------%
            dv = pds.spikeserver.connect(dv);
            
            % From help PsychDataPixx:
            % Timestamping is disabled by default (mode == 0), as it incurs a bit of
            % computational overhead to acquire and log timestamps, typically up to 2-3
            % msecs of extra time per 'Flip' command.
            % Buffer is collected at the end of the expeiment!
            PsychDataPixx('LogOnsetTimestamps', 2);%2
            PsychDataPixx('ClearTimestampLog');
            
    
            % Initialize Datapixx for Dual CLUTS
            dv = pds.datapixx.init(dv);
            
            pds.keyboard.setup();
    

    %% Last chance to check variables
    if(dv.trial.pldaps.pause.type==1 && dv.trial.pldaps.pause.preExperiment==true) %0=don't,1 is debugger, 2=pause loop
        dv  %#ok<NOPRT>
        disp('Ready to begin trials. Type return to start first trial...')
        keyboard %#ok<MCKBD>
    end
 
    %%%%start recoding on all controlled components this in not currently done here
    % save timing info from all controlled components (datapixx, eyelink, this pc)
    dv = beginExperiment(dv);

    % disable keyboard
    ListenChar(2)
    HideCursor
    
    %save defaultParameters as trial 0
    trialNr=0;
    dv.trial.pldaps.iTrial=0;
    dv.trial=mergeToSingleStruct(dv.defaultParameters);
    result = saveTempFile(dv); 
    if ~isempty(result)
        disp(result.message)
    end
    
    
    %now setup everything for the first trial
   
%     dv.defaultParameters.pldaps.iTrial=trialNr;
    
    %we'll have a trialNr counter that the trial function can tamper with?
    %do we need to lock the defaultParameters to prevent tampering there?
    levelsPreTrials=dv.defaultParameters.getAllLevels();
%     dv.defaultParameters.addLevels(dv.conditions(trialNr), {['Trial' num2str(trialNr) 'Parameters']});
    
    %for now all structs will be in the parameters class, first
    %levelsPreTrials, then we'll add the condition struct before each trial.
%     dv.defaultParameters.setLevels([levelsPreTrials length(levelsPreTrials)+trialNr])
%     dv.defaultParameters.pldaps.iTrial=trialNr;
%     dv.trial=mergeToSingleStruct(dv.defaultParameters);
    
    %only use dv.trial from here on!
    
    %% main trial loop %%
    while dv.trial.pldaps.iTrial < dv.trial.pldaps.finish && dv.trial.pldaps.quit~=2
        
        if dv.trial.pldaps.quit == 0
            
           %load parameters for next trial and lock defaultsParameters
           trialNr=trialNr+1;
           dv.defaultParameters.addLevels(dv.conditions(trialNr), {['Trial' num2str(trialNr) 'Parameters']});
           dv.defaultParameters.setLevels([levelsPreTrials length(levelsPreTrials)+trialNr]);
           dv.defaultParameters.pldaps.iTrial=trialNr;
           dv.trial=mergeToSingleStruct(dv.defaultParameters);
           dv.defaultParameters.setLock(true);
            
           % run trial
           dv = feval(dv.trial.pldaps.trialMasterFunction,  dv);
            
           %unlock the defaultParameters
           dv.defaultParameters.setLock(false); 
            
           %save tmp data
           result = saveTempFile(dv); 
           if ~isempty(result)
               disp(result.message)
           end
                      
           %store the difference of the trial struct to .data
           dTrialStruct=getDifferenceFromStruct(dv.defaultParameters,dv.trial);
           dv.data{trialNr}=dTrialStruct;
           
           %advance to next trial
%            if(dv.trial.pldaps.iTrial ~= dv.trial.pldaps.finish)
%                 %now we add this and the next Trials condition parameters
%                 dv.defaultParameters.addLevels(dv.conditions(trialNr), {['Trial' num2str(trialNr) 'Parameters']},[levelsPreTrials length(levelsPreTrials)+trialNr]);
%                 dv.defaultParameters.pldaps.iTrial=trialNr;
%                 dv.trial=mergeToSingleStruct(dv.defaultParameters);
%            else
%                 dv.trial.pldaps.iTrial=trialNr;
%            end
%            
%            if isfield(dTrialStruct,'pldaps')
%                if isfield(dTrialStruct.pldaps,'finish') 
%                     dv.trial.pldaps.finish=dTrialStruct.pldaps.finish;
%                end
%                if isfield(dTrialStruct.pldaps,'quit') 
%                     dv.trial.pldaps.quit=dTrialStruct.pldaps.quit;
%                end
%            end
            
        else %dbquit ==1 is meant to be pause. should we halt eyelink, datapixx, etc?
            %create a new level to store all changes in, 
            %load only non trial paraeters
            pause=dv.trial.pldaps.pause.type;
            dv.trial=dv.defaultParameters;
            
            dv.defaultParameters.addLevels({struct}, {['PauseAfterTrial' num2str(trialNr) 'Parameters']});
            dv.defaultParameters.setLevels([levelsPreTrials length(dv.defaultParameters.getAllLevels())]);
            
            if pause==1 %0=don't,1 is debugger, 2=pause loop
                ListenChar(0);
                ShowCursor;
                dv.trial
                disp('Ready to begin trials. Type return to start first trial...')
                keyboard %#ok<MCKBD>
                dv.trial.pldaps.quit = 0;
                ListenChar(2);
                HideCursor;
            elseif pause==2
                pauseLoop(dv);
            end           
%             pds.datapixx.refresh(dv);

            %now I'm assuming that nobody created new levels,
            %but I guess when you know how to do that
            %you should also now how to not skrew things up
            allStructs=dv.defaultParameters.getAllStructs();
            if(~isequal(struct,allStructs{end}))
                levelsPreTrials=[levelsPreTrials length(allStructs)]; %#ok<AGROW>
            end
        end
        
    end
    
    %make the session parameterStruct active
    dv.defaultParameters.setLevels(levelsPreTrials);
    dv.trial = dv.defaultParameters;
    
    % return cursor and command-line control
    ShowCursor;
    ListenChar(0);
    Priority(0);
    
    dv = pds.eyelink.finish(dv);
    dv = pds.spikeserver.disconnect(dv);
    if(dv.defaultParameters.datapixx.use)
        dv.defaultParameters.datapixx.timestamplog = PsychDataPixx('GetTimestampLog', 1);
    end
    
    
    if ~dv.defaultParameters.pldaps.nosave
        [structs,structNames] = dv.defaultParameters.getAllStructs();
        
        PDS=struct;
        PDS.initialParameters=structs(levelsPreTrials);
        PDS.initialParameterNames=structNames(levelsPreTrials);
        PDS.initialParametersMerged=mergeToSingleStruct(dv.defaultParameters); %too redundant?
        PDS.conditions=structs((max(levelsPreTrials)+1):end);
        PDS.conditionNames=structNames((max(levelsPreTrials)+1):end);
        PDS.data=dv.data; 
        save(fullfile(dv.defaultParameters.session.dir, dv.defaultParameters.session.file),'PDS','-mat')
    end
    
    
    Screen('CloseAll');
    sca;
    
    
catch me
    sca
    
    % return cursor and command-line control
    ShowCursor
    ListenChar(0)
    disp(me.message)
    
    nErr = size(me.stack); 
    for iErr = 1:nErr
        fprintf('errors in %s line %d\r', me.stack(iErr).name, me.stack(iErr).line)
    end
    fprintf('\r\r')
    keyboard    
end

end
%we are pausing, will create a new defaultParaneters Level where changes
%would go.
function pauseLoop(dv)
        ShowCursor;
        ListenChar(1);
        while(true)
            %the keyboard chechking we only capture ctrl+alt key presses.
            [dv.trial.keyboard.pressedQ,  dv.trial.keyboard.firstPressQ]=KbQueueCheck(); % fast
            if dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.Lctrl)&&dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.Lalt)
                %D: Debugger
                if  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.dKey) 
                    disp('stepped into debugger. Type return to start first trial...')
                    keyboard %#ok<MCKBD>

                %E: Eyetracker Setup
                elseif  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.eKey)
                    try
                       if(dv.trial.eyelink.use) 
                           pds.eyelink.calibrate(dv);
                       end
                    catch ME
                        display(ME);
                    end

                %M: Manual reward
                elseif  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.mKey)
                    if dv.trial.datapixx.use
                        pds.datapixx.analogOut(dv.trial.stimulus.rewardTime)
                        pds.datapixx.flipBit(dv.trial.event.REWARD);
                    end
                    dv.trial.ttime = GetSecs - dv.trial.trstart;
                    dv.trial.stimulus.timeReward(:,dv.trial.iReward) = [dv.trial.ttime dv.trial.stimulus.rewardTime];
                    dv.trial.stimulus.iReward = dv.trial.iReward + 1;
                    PsychPortAudio('Start', dv.trial.sound.reward);

                %P: PAUSE (end the pause) 
                elseif  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.pKey)
                    dv.trial.pldaps.quit = 0;
                    ListenChar(2);
                    HideCursor;
                    break;

                %Q: QUIT
                elseif  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.qKey)
                    dv.trial.pldaps.quit = 2;
                    break;
                
                %X: Execute text selected in Matlab editor
                elseif  dv.trial.keyboard.firstPressQ(dv.trial.keyboard.codes.xKey)
                    activeEditor=matlab.desktop.editor.getActive; 
                    if isempty(activeEditor)
                        display('No Matlab editor open -> Nothing to execute');
                    else
                        if isempty(activeEditor.SelectedText)
                            display('Nothing selected in the active editor Widnow -> Nothing to execute');
                        else
                            try
                                eval(activeEditor.SelectedText)
                            catch ME
                                display(ME);
                            end
                        end
                    end
                    
                    
                end %IF CTRL+ALT PRESSED
            end
            pause(0.1);
        end

end
