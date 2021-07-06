function p = pmTemplate_matrixModuleGabor(p, state, sn)
% function p = <yourPackage>.pmTemplate_matrixModuleGabor(p, state, sn)
% 
% A "matrixModule" --short for 'condition matrix module'-- is a special type
% of PLDAPS module used by the PLDAPS condMatrix to apply & update parameters
% that define the condition(s) presented on every trial
% - the condMatrix identifies all matrixModules by the flag:
%       p.trial.(sn).stateFunction.matrixModule == 1
% 
% Template for a "matrixModule"using modular PLDAPS experimental design
% - "pm" in filename is convention for PLDAPS module
%   - *ALL* PLDAPS modules must have these--and only these--3 inputs:
%     (p, state, sn)
% 
% INPUTS:
%   [p]      PLDAPS object
%   [state]  Numeric value of the current PLDAPS state, e.g. trialSetup, frameUpdate, frameDraw, etc. (see p.trial.pldaps.trialStates)
%   [sn]     String name of the module field within p.trial struct, i.e. p.trial.(sn);
% 
% 
% --- USAGE ---
% Initialization & formatting of the [sn] input is a little distinct from a standard module.
% Because matrixModules are actually a family of 1 or more modules that will
% be automatically updated on each trial, its helpful to have a "base module"
% field that can store certain parameters that are shared by all matrixModules.
% To facilitate this, all matrixModules should be created with a 2 digit index
% number appended to the end of the base module name.
% e.g.  'gabor' matrix modules might share a common size, spatial frequency, origin
%       location, etc
%       - place these common parameters inside p.trial.gabor
%       - create the desired number of matrix modules that comprise every trial
%         as p.trial.gabor01, p.trial.gabor02...p.trial.gabor##
%       - the code of the matrix module can then extract the base name [snBase]
%         by simply trimming off the last two characters of the [sn] input
%
% --- inside  yourPackage.doMyExperiment.m
% 
%     %% Initialization
%     % Initialize matrixModule(s) inside your "do" experiment wrapper by first setting up
%     % [sn] as a normal struct field of your PLDAPS settings struct [pss]:
% 
%     %% Setup shared matrixModule parameters & components
%     sn = 'gabors';
% 
%     % Base params field (non-module) stores params shared across all matrixModule instances
%     pss.(sn).stimCtr = [0, 0, 0]; % global stimulus center location; visual degrees relative to screen center
%     pss.(sn).gabTf = 3; % gabor temporal frequency
%     pss.(sn).gabFwhm = 2; % full-width half-max size of each gabor (visual degrees)
% 
%     pss.(sn).centerOnScreen = true; % flag to use procedural gabor creation that has it's origin at screen center (rather than upper left corner)
%     pss.(sn).trackMouse = 1;        
% 
%     % Create a temporary module that will be used to create each of your matrixModules
%     tmpModule = pldapsModule('modName',sn, 'name','templates.pmTemplate_matrixModuleGabor', 'matrixModule',true, 'order',10); 
% 
%     % Apply any adjustments or defaults necessary in the matrixModule
%     tmpModule.use = true;   % whether or not PLDAPS state machine should use the module when iterating through experiment states
%                             % - important during basic initializaiton states like .experimentPreOpenScreen, 
%     tmpModule.on = false;   % whether or not the module should be seen by your behavioral state machine (see modularDemo.pmBase.m)
%   
%     % Populate shared stimulus parameters 
%     tmpModule.gabContrast = 1;
%     tmpModule.ngabors = 1;
%     % Some parameters are direclty copied over from the stimulus struct:
%     % - NOTE that its just .(sn) during setup, but within the module code this shared location will become .(snBase)
%     tmpModule.centerOnScreen = pss.(sn).centerOnScreen;
%     tmpModule.gabTf = pss.(sn).gabTf; % drift rate (Hz)
%     tmpModule.gabFwhm = pss.(sn).gabFwhm;
%     tmpModule.gabSf = 1.3/tmpModule.gabFwhm;    % V1 median bandwidth= ~1.4 (== 1.0/fwhm) --De Valois 1982
% 
%     % Initialize parameters to be updated by condMatrix prior to each trial
%     tmpModule.pos = zeros(2, tmpModule.ngabors); % [x,y] position of stimulus; visual degrees **relative to .stimCenter** (see code below)
%     tmpModule.dir = zeros(1, tmpModule.ngabors);
%     
%     % Select module stimulus type
%     % - Not necessary, but sometimes useful to have multiple 'flavors' of a particular module
%     % - Makes experiment code more robust/universal
%     % - Biggest benefit is ability to use common analysis code across experiments with the same general stimulus elements (i.e. gabors, dots, etc)
%     tmpModule.type = 'steerX';    %'polarTrack';
%     
%     % Stimulus onset timing & n-reps per trial
%     ncopies = 1;
%     tmpModule.isi = .0; % blank interstimulus interval between each matrix module
%     % stimDur==total "stimulus on" duration within the trial
%     % - should be setup while initializing your behavioral module (e.g. modularDemo.pmBase.m)
%     stimModuleDur = stimDur/ncopies;
% 
% 
%     %% Create matrixModules
%     % Make duplicate/indexed stim modules for each repetition w/in a trial
%     matrixModNames = {};
%     for i = 1:ncopies
%         mN = sprintf('%s%0.2d', sn, i); % indexed module name
%         pss.(mN) = tmpModule;   % apply the temporary module
%         pss.(mN).stateFunction.modName = mN; % Update it's internalized name field
%         % pss.(mN).stateFunction.matrixModule = true; % ...we already set this flag using a param-value pair input to pldapsModule()
%         
%         % timing: each module [onset, offset] time in sec (relative to STIMULUS state start)
%         basedur = (i-1)*stimModuleDur;
%         pss.(mN).modOnDur = [0, stimModuleDur-tmpModule.isi] +basedur; % subtract isi from offset
% 
%         % Compile module names [for use later within experiment wrapper code (...if needed)]s
%         matrixModNames{i} = mN;
%     end
% 
% ---
% 
% Dependencies:
%     visBasics.     .m
%     pds.applyDefaults.m
%     PLDAPS/SupportFunctions/Utils/gaborShader.*          
%
% 2018-xx-xx  TBC  Wrote it
% 2020-10-16  TBC  Gabor positioning update:
%                  - NEW [.centerOnScreen] flag to use screen center as origin,
%                  -- destRect texture coords now consistent with Screen('DrawDots'..)
% 2021-03-10  TBC  Created & commented for basic template use
% ---------------------------------------------------------------- 
% Stimulus position & screen center:
% Stimulus rendering coordinates of PTB & OpenGL can be confusing/conflicting from
% time to time.
% PTB's default places the origin in the upper left corner of the screen, with
% the positive Y-axis in the downward direction.
% 
% OpenGL (*logically*) places the origin at the center of the screen, with
% the positive Y-axis in the upward direction.
% (+x,+y) is in the upper right quadrant
% (-x,-y) is in the lower left quadrant.
% 
% This module attempts to rectify this discrepancy by aligning to the OpenGL
% coordinate frame. Because centering on screen involves both a shift & y-axis
% flip (relative to PTB pixel-centric coordinates) this leads to some crufty
% confusing code, but we've tried to minimize & streamline it as much as possible.
%   -- T.Czuba, Nov. 2020
% 
% ----------------------------------------------------------------

snBase = sn(1:end-2);
% base name by stripping off matrix module number (2 digits)
%  ...my own dumb coding. --TBC

% stimulus relevant PLDAPS modules
snBehav = p.trial.pldaps.modNames.behavior{1}; % behavioral control module
snFix = p.trial.pldaps.modNames.currentFix{1}; % current fixation module

% Flag operations that should only be done once across ALL matrixModule instances
oneTimers = contains(sn(end-1:end),'01');

switch state
    % FRAME STATES
    case p.trial.pldaps.trialStates.frameUpdate
            % [T]rack mouse position
            if  p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.tKey)
                p.trial.(snBase).trackMouse = true;
                
            elseif p.trial.keyboard.firstPressQ(p.trial.keyboard.codes.rKey)
                % [R]emain in the same position (..?)
                p.trial.(snBase).trackMouse = false;
                try
                    mousePos = [1;-1].*(p.trial.mouse.cursorSamples(1:2,p.trial.mouse.samples)-p.trial.display.ctr(1:2)');
                    fprintf('Cursor Pos: %3.3f, %3.3f deg\n', pds.px2deg(mousePos, p.trial.display.viewdist, p.trial.display.px2w));
                end
                
            end
            
            if p.trial.(sn).on
                if p.trial.(snBase).trackMouse
                    % dynamically update stimulus position based on [change in] mouse position
                    updateGaborPosition(p, sn);
                end
                
                % Check for trial completion:
                % - compute stimulus position relative to fixation target
                % - if w/in desired range, flag as response/complete
                % - here, we utilize existing code for checking eye position relative to fixation point (snFix)
                %   but instead of letting it poll the current .eyeX & .eyeY positions,
                %   we pass it the current stimulus position (established w/in updateGaborPosition) as input.
                p.trial.(snBehav).success = pldaps.checkFixation(p, snFix, p.trial.(sn).pos(1:2)');
            end

            
    case p.trial.pldaps.trialStates.framePrepareDrawing
        if p.trial.(sn).on
            % Record stim position for this frame:   [xyz, element, frame]
            p.trial.(sn).posFrames(:,:,p.trial.iFrame) = p.trial.(sn).pos;

            % update phase with shared TF parameter
            phaseStep = p.trial.(snBase).gabTf/p.trial.display.frate*360;
                        
            % drift gabor phase w/in procedural drawing parameters
            for i = p.trial.display.bufferIdx
                % SAME direction in either eye
                p.trial.(sn).Gpars(1,:,i+1) = p.trial.(sn).Gpars(1,:,i+1) + phaseStep;
                % % opposite direction in either eye would be:
                % p.trial.(sn).Gpars(1,:,i+1) = p.trial.(sn).Gpars(1,:,i+1) + phaseStep * (-i*2+1);
            end
                        
        end
        
        
    case p.trial.pldaps.trialStates.frameDraw
        if p.trial.(sn).on
            % must change blend function for procedural gabor rendering
            % ...only one BlendFunction to set, not separate for each stereobuffer.
            [prevSrcRect, prevDestRect] = Screen('BlendFunction', p.trial.display.ptr, 'GL_ONE', 'GL_ONE');
            drawTheGabors(p, sn);
            % reset previous blend functions
            Screen('BlendFunction', p.trial.display.ptr, prevSrcRect, prevDestRect);
        end
        
        
	% TRIAL STATES        
    case p.trial.pldaps.trialStates.trialPrepare
        % only do here what can't be done in conditions matrix
        trialPrepare(p, sn);
        
        
    % EXPT STATES
    case p.trial.pldaps.trialStates.experimentPreOpenScreen
        initParams(p, sn);

        
    case p.trial.pldaps.trialStates.experimentPostOpenScreen        
        % Create procedural textures
        %   -- input params assumed to be in visual degrees
        p = modularDemo.makeGaborsCentered(p, sn);
        
        
    case p.trial.pldaps.trialStates.experimentCleanUp
        Screen('Close', p.trial.(sn).gabTex);
        
end


%% Nested functions
    
    function p = updateGaborPosition(p, sn)
                
        if p.trial.(snBase).trackMouse % && p.trial.mouse.samples>0
            % index of two most recent mouse position samples
            ii = p.trial.mouse.samples + [-1,0];
            ii(ii<1) = 1; % must be valid indices
            % compute difference (invert Y)
            mouseDelta = diff(p.trial.mouse.cursorSamples(:,ii)' .* [1,-1])';
            % & convert to visual degrees
            mouseDelta = mousePx2Deg(p, mouseDelta, 0);
            % accumulate this [xy] offset throughout trial
            % - this param must be w/in each module, else the accumulated position offset
            %   will be efffectively multiplied by n-matrixModules(!)
            p.trial.(sn).stimOffset(1:2) = p.trial.(sn).stimOffset(1:2) + mouseDelta(1:2);
            
            stimCtr = p.trial.(snBase).stimCtr + p.trial.(sn).stimOffset;
            
        else
            % use .eyeDelta?
            try
                % - this value [may be] precomputed from whatever is being tracked by the .tracking module
                p.trial.(sn).stimOffset(1:2) = p.trial.(sn).stimOffset(1:2) + p.trial.eyeDelta(1:2);
            end
        end
        % stimCtr = p.trial.(snBase).stimCtr;
        stimCtr = p.trial.(snBase).stimCtr + p.trial.(sn).stimOffset;
        
        
        % Calculate/update gabor positions of THIS module (in deg)
        for i = 1:2 % [x,y]
            p.trial.(sn).pos(i,:) = (p.trial.(sn).stimPos(i) + stimCtr(i));
            %NOTE: mouse samples don't exist yet in PLDAPS [trialPrepare] state (...not prior to sampling in pldapsDefaultTrial>>frameUpdate state)
        end
        
        % Convert pos visual degrees to pixels
        p.trial.(sn).pos = p.trial.(sn).pos  .* p.trial.display.ppd;
        % NOTE on deg2pixel conversion:
        % Non-projection method using simple pixel-per-degree(.ppd) approximation
        % Alternatively, could use projective geometry with:
        %       p.trial.(sn).pos = pds.px2deg(p.trial.(sn).pos, p.trial.display.viewdist, p.trial.display.px2w);
        % BUT must take extra care that projective conversion is **only applied once**
        
        % Gabor bounding rect:
        gabRect = p.trial.(sn).texRect;
        
        % update gabor rects to this positon
        p.trial.(sn).gabRects = CenterRectOnPointd( gabRect, p.trial.(sn).pos(1,:)', p.trial.(sn).pos(2,:)');

    end

%% trialPrepare
    function trialPrepare(p, sn)

        % Update gabor parameters & conversions (covers params set in procedural gabor creation function:  glDraw.makeGabors.m)
        % - necessary for gabor shape/size/sf parameter flexibility, and viewing distance dependence
        gabSd = p.trial.(snBase).gabFwhm ./ sqrt(8*log(2));
        if ~isequal(p.trial.(sn).gabSd, gabSd)
            p.trial.(sn).gabSd = gabSd;
        end
        p.trial.(sn).gabPixCycle    = p.trial.(sn).gabSf / p.trial.display.ppd;   % pix/deg * deg/cycle = pix/cycle
        p.trial.(sn).gabPixels      = 8 * p.trial.(sn).gabSd .* p.trial.display.ppd;
        p.trial.(sn).texSz          = ceil(p.trial.(sn).gabPixels);
        p.trial.(sn).texRect        = [0 0 1 1] * p.trial.(sn).texSz(:);
        
        p.trial.(sn).Gpars(2,:)     = p.trial.(sn).gabPixCycle;
        p.trial.(sn).Gpars(3,:)     = p.trial.(sn).gabSd * p.trial.display.ppd; % std of gaussian hull (in pixels)
                
        % update gabor rects based on position
        p = updateGaborPosition(p, sn);

        % Record position on a frame-by-frame basis:
        % - for analysis of data collected during dynamic stim update (.trackMouse == true)
        p.trial.(sn).posFrames = nan( [size(p.trial.(sn).pos), ceil(p.trial.pldaps.maxFrames)] );
        
        
        % Randomize initial phase
        p.trial.(sn).Gpars(1,:) = rand([1, p.trial.(sn).ngabors*numel(p.trial.display.bufferIdx)])*360;
        
        % If stereo, establish interocular phase
        if numel(p.trial.display.bufferIdx)>1
%             switch p.trial.(sn).binoPhaseMode
%                 case 0
                    % ZERO interocular phase difference (i.e. in plane of fixation)
                    p.trial.(sn).Gpars(1,:,2) = p.trial.(sn).Gpars(1,:,1);
%                 case 1
%                     % MATCHED random initial phase offset for all gabors
%                     ioPhaseDiff = diff(squeeze(p.trial.(sn).Gpars(1,1,:))); % use difference of first pair (...gives 1:1 link to fully random version)
%                     p.trial.(sn).Gpars(1,:,2) = p.trial.(sn).Gpars(1,:,1) + ioPhaseDiff;
%                 case 2
%                     % FULLY RANDOM initial phase offsets
%                     % do nothing
%             end
        end
        
    end % end trialPrepare


    %% initParams
    % Initialize default module parameters
    function initParams(p, sn)
        % list of default parameters
        def = struct(...
            'on', true ...
            ,'gabSd', .43 ...
            ,'gabSf', 2 ...
            ,'gabContrast', 1 ...
            ,'ngabors', 1 ...
            ,'pos', [0 0 0]' ...
            ,'dir', 0 ...
            ,'centerOnScreen', true ...
            ,'stimOffset', [0 0 0]' ...
            ,'markerRects', [] ...
            ,'drawMarkers', false ...
            ,'trackMouse', false ...
            );
        
        % "pldaps requestedStates" for this module
        % - done here so that modifications only needbe made once, not everywhere this module is used
        rsNames = {'frameUpdate', 'framePrepareDrawing', 'frameDraw', ...
                   'trialPrepare', ...'trialCleanUpandSave', ...
                   'experimentPreOpenScreen', 'experimentPostOpenScreen', 'experimentCleanUp'};
               
        % Apply default params
        p.trial.(sn) = pds.applyDefaults(p.trial.(sn), def);
        p = pldapsModuleSetStates(p, sn, rsNames);
        
    end % end initParams
    

    %% drawTheGabors
    function drawTheGabors(p, sn)
        % .texRectCtr is centered on screen .display.
        if p.trial.(sn).on
            for i = p.trial.display.bufferIdx
                Screen('SelectStereoDrawBuffer',p.trial.display.ptr, i);
                Screen('DrawTextures', p.trial.display.ptr, p.trial.(sn).gabTex, [], p.trial.(sn).gabRects', p.trial.(sn).dir(:,min([i+1,end])), [], [], [], [], kPsychDontDoRotation, p.trial.(sn).Gpars(:,:,i+1));
            end
        end        
    end % end drawTheGabors
    

%% mousePx2Deg
    function degOut = mousePx2Deg(p, pxIn, centerPxIn)            
        if nargin<2
            % get most recent mouse pos
            pxIn = p.trial.mouse.cursorSamples(:,p.trial.mouse.samples);
        end
        if nargin<3 || centerPxIn
            % Default: correct mouse pos to upright screen-centered coordinates (i.e. glDraw coordinates)
            pxIn = [1;-1] .* (pxIn - p.trial.display.ctr(1:2)');
        end
        
        % simple conversion using pixel-per-degree(.ppd) approximation
        degOut = pxIn ./ p.trial.display.ppd;
        
    end %mousePx2Deg


end