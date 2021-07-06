function p =  init(p)
%pds.datapixx.init    initialize Datapixx at the beginning of an experiment.
% p =  pds.datapixx.init(p)
%
% pds.datapixx.init is a function that intializes the DATAPIXX, preparing it for
% experiments. Critically, the PSYCHIMAGING calls sets up the dual CLUTS
% (Color Look Up Table) for two screens.  These two CLUTS are in the
% p.trial.display substruct
% INPUTS
%       p.trial [struct] - main pldaps display variables structure
%           .dispplay [struct] - required display
%               .ptr         - pointer to open PTB window
%               .useOverlay  - (0,1,2) for whether to use CLUT overlay
%                       0 - no overlay
%                       1 - Datapixx overlay (if datapixx.use is off, draws
%                                             overlay version)
%                       2 - Software overlay (requires stereomode setup)
%               .gamma.table - required ( can be a linear gamma table )
%               .humanCLUT   [256 x 3] human color look up table
%               .monkeyCLUT  [256 x 3] monkey color look up table
% OUTPUTS
%       p [modified]
%           .trial.display.overlayptr - pointer to overlay window added
%
%           if useOverlay == 1, pds.datapixx.init opens an overlay pointer
%           with p.trial.display.monkeyCLUT as the color look up table for one
%           datapixx monitor and p.trial.display.humanCLUT for the other monitor.
%
% Datapixx is Opened and set to default settings for PLDAPS


% 2011       kme wrote it
% 12/03/2013 jly reboot for version 3
% 2014       adapt to version 4.1
% 2016       jly add software overlay
% 2018-05-02 tbc  CombinedClut now always applied via Datapixx mex, instead of via Screen

global dpx GL;

%% initialize Datapixx
if p.trial.datapixx.use
    
    if ~Datapixx('IsReady')
        Datapixx('Open');
    end
    
    % From help PsychDataPixx:
    % Timestamping is disabled by default (mode == 0), as it incurs a bit of
    % computational overhead to acquire and log timestamps, typically up to 2-3
    % msecs of extra time per 'Flip' command.
    % Buffer is collected at the end of the expeiment!
    if p.trial.datapixx.LogOnsetTimestampLevel==2
        % .LogOnsetTimestampLevel==2 is NOT RECOMMENDED (very slow).
        % When ==2, ALL Screen('Flip') times will be logged w/in the datapixx box,
        % then transferred at the end of the experiment
        % WHEN ==1 (PLDAPS default), a high-precision Datapixx timestamp will only be
        % created ONCE at the start of each trial (w/in pldapsDefaultTrial.m)
        % In either case, Datapixx timestamp log is transferred at the completion of an
        % experiment session (w/in pldaps.run.m), and appear in:
        %   PDS.baseParams.datapixx.timestamplog
        PsychDataPixx('LogOnsetTimestamps',p.trial.datapixx.LogOnsetTimestampLevel);
        warnStr = sprintf([fprintLineBreak('!'), ...
                           '~!~\t Warning: Your PLDAPS .datapixx setup is configured for\n', ...
                           '~!~\t additional high-precision timestamp LOGGING OF EVERY FRAME\n', ...
                           '~!~\t\t .datapixx.LogOnsetTimestampLevel == 2;\n', ...
                           '~!~\t Acquisition & transfer of this tends to be a very time-intensive \n', ...
                           '~!~\t (& occasionally error prone) process, thus is GENERALLY NOT RECOMMENDED.\n', ...
                           '~!~\t\n~!~\t See help PsychDataPixx  (section on ''LogOnsetTimestamps'')\n', ...
                           fprintLineBreak('!'), ...
                           'Default behavior is to only acquire one high-precision timestamp at the start of each trial,\n', ...
                           'that can later be used to detect any stimulus onset discrepancies in the standard PTB timestamps\n', ...
                           'Default:\n', ...
                           '\t .datapixx.LogOnsetTimestampLevel == 1;\n', ...
                           '\t One high-precision timestamp/trial will appear (at end of experiment) in:\n', ...
                           '\t\t PDS.baseParams.datapixx.timestamplog.\n', ...
                           fprintLineBreak('-'), ...
                           'As always, standard PTB timestamps for every frame appear (as they happen) in:\n', ...
                           '\t\t p.data{}.timing.flipTimes\n', ...
                           fprintLineBreak('-')]);
        fprintf(2, warnStr)
    end
    PsychDataPixx('ClearTimestampLog');
    
    %set getPreciseTime options, see testsuite/pldapsTimingTests for
    %details
    if isfield(p.trial.datapixx.GetPreciseTime, 'syncmode') && ~isempty(p.trial.datapixx.GetPreciseTime.syncmode)
        dpx.syncmode=p.trial.datapixx.GetPreciseTime.syncmode; %1,2,3
    end
    if ~isempty(p.trial.datapixx.GetPreciseTime.maxDuration)
        dpx.maxDuration=p.trial.datapixx.GetPreciseTime.maxDuration;
    end
    if ~isempty(p.trial.datapixx.GetPreciseTime.optMinwinThreshold)
        dpx.optMinwinThreshold=p.trial.datapixx.GetPreciseTime.optMinwinThreshold;
    end
    
    if Datapixx('IsPropixx')
        %this might not work reliably
        if ~Datapixx('IsPropixxAwake')
            Datapixx('SetPropixxAwake');
        end
        Datapixx('EnablePropixxLampLed');
        % Datapixx('RegWrRd');
        
        if p.trial.datapixx.enablePropixxRearProjection
            Datapixx('EnablePropixxRearProjection');
        else
            Datapixx('DisablePropixxRearProjection');
        end
        % Datapixx('RegWrRd');
        
        if p.trial.datapixx.enablePropixxCeilingMount
            Datapixx('EnablePropixxCeilingMount');
        else
            Datapixx('DisablePropixxCeilingMount');
        end
        % Datapixx('RegWrRd');
        
        if isfield(p.trial.datapixx,'rb3d') && p.trial.datapixx.rb3d
            if p.trial.display.useOverlay==1
                Datapixx('SetVideoMode', 9); % Enable the overlay on the PPX CTRL
            end
            Datapixx('SetPropixxDlpSequenceProgram', 1); %, 1);
            Datapixx('RegWrRd');
        end
        % Datapixx('RegWrRd');

    end
    
    p.trial.datapixx.info.DatapixxFirmwareRevision = Datapixx('GetFirmwareRev');
    p.trial.datapixx.info.DatapixxRamSize = Datapixx('GetRamSize');
       
    %%% Open Datapixx and get ready for data aquisition %%%
    Datapixx('StopAllSchedules');
    Datapixx('DisableDinDebounce');
    Datapixx('EnableAdcFreeRunning');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('SetDoutValues',0);
    Datapixx('RegWrRd');
    
    %start adc data collection if requested
    pds.datapixx.adc.start(p);
    
    % Create field for pending strobed word values
    p.trial.datapixx.strobeQ = [];
    
end

if p.trial.display.useOverlay==1 % Datapixx overlay
    if p.trial.datapixx.use
        %disp('****************************************************************')
        fprintLineBreak
        disp('Initializing Datapixx Overlay')
        %disp('Combining color look up tables that can be found in')
        %disp('dv.disp.humanCLUT and dv.disp.monkeyCLUT')
        %disp('****************************************************************')
        
        % Set overlay transparency color to background, ensuring it matches post-gamma correction
        bgColor=p.trial.display.bgColor;
        if isField(p.trial, 'display.gamma.table')
            bgColor = interp1(linspace(0,1,256),p.trial.display.gamma.table(:,1), p.trial.display.bgColor);
        elseif isField(p.trial, 'display.gamma.power')
            % outcolor = incolor ^ EncodingGamma.
            bgColor =  p.trial.display.bgColor .^ p.trial.display.gamma.power;
        end
        Datapixx('SetVideoClutTransparencyColor', bgColor);
        Datapixx('EnableVideoClutTransparencyColorMode');
        Datapixx('RegWrRd');
        
        if p.trial.display.switchOverlayCLUTs
            combinedClut = [p.trial.display.humanCLUT; p.trial.display.monkeyCLUT];
        else
            combinedClut = [p.trial.display.monkeyCLUT; p.trial.display.humanCLUT];
        end
        %%% Gamma correction for dual CLUT %%%
        % check if gamma correction has been run on the window pointer
        if isField(p.trial, 'display.gamma.table')
            % get size of the combiend CLUT. It should be 512 x 3 (two 256 x 3 CLUTS
            % on top of eachother).
            sc = size(combinedClut);
            
            % use sc to make a vector of 8-bit color steps from 0-1
            x = linspace(0,1,sc(1)/2);
            % use the gamma table to lookup what the values should be
            y = interp1(x,p.trial.display.gamma.table(:,1), combinedClut(:));
            % reshape the combined clut back to 512 x 3
            combinedClut = reshape(y, sc);
        elseif isField(p.trial, 'display.gamma.power')            
            combinedClut=combinedClut .^ p.trial.display.gamma.power;
            
        end

        % Retrieve/update extended Datapixx settings ("VideoStatus")
        p.trial.datapixx.videoStatus = Datapixx('GetVideoStatus');
        
        % Apply the clut directly to the Datapixx
        Datapixx('SetVideoClut', combinedClut)
        Datapixx('RegWrRd');
        % formerly done through Screen, but multiple systems now require the clut to be set by
        % Datapixx mex directly for TransparencyColors to work (e.g. ProPixx & ViewPixx3D displays as of 2018)

    end
    
elseif p.trial.display.useOverlay==2 % software overlay
    fprintLineBreak
    disp('Initializing PLDAPS Software Overlay')
    %assign transparency color
    bgColor=p.trial.display.bgColor;
    glUniform3f(glGetUniformLocation(p.trial.display.shader, 'transparencycolor'), bgColor(1), bgColor(2), bgColor(3));
    
    if p.trial.display.switchOverlayCLUTs
        combinedClut = [p.trial.display.humanCLUT; p.trial.display.monkeyCLUT];
    else
        combinedClut = [p.trial.display.monkeyCLUT; p.trial.display.humanCLUT];
    end

    % assign values to look up textures
    % Bind relevant texture object:
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, p.trial.display.lookupstexs(1));
    % Set filters properly: Want nearest neighbour filtering, ie., no filtering
    % at all. We'll do our own linear filtering in the ICM shader. This way
    % we can provide accelerated linear interpolation on all GPU's with all
    % texture formats, even if GPU's are old:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST)
    % Want clamp-to-edge behaviour to saturate at minimum and maximum
    % intensity value, and to make sure that a pure-luminance 1 row clut is
    % properly "replicated" to all three color channels in rgb modes:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);
    % Assign lookuptable data to texture:
    n=size(p.trial.display.humanCLUT, 1);
    m=size(p.trial.display.humanCLUT, 2);
    glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, p.trial.display.internalFormat,n,m, 0,GL.LUMINANCE, GL.FLOAT, single(combinedClut(1:n,:)));
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);

    %#2
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, p.trial.display.lookupstexs(2));
    % Set filters properly: Want nearest neighbour filtering, ie., no filtering
    % at all. We'll do our own linear filtering in the ICM shader. This way
    % we can provide accelerated linear interpolation on all GPU's with all
    % texture formats, even if GPU's are old:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST)
    % Want clamp-to-edge behaviour to saturate at minimum and maximum
    % intensity value, and to make sure that a pure-luminance 1 row clut is
    % properly "replicated" to all three color channels in rgb modes:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);
    % Assign lookuptable data to texture:
    glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, p.trial.display.internalFormat, n, m, 0, GL.LUMINANCE, GL.FLOAT, single(combinedClut(n+1:end,:)));
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
end


% 3D crosstalk correction when using Propixx Rb3d  (...can operate independent of overlay)
if isfield(p.trial.datapixx, 'crosstalk') && any(p.trial.datapixx.crosstalk(:))
    Datapixx('SetPropixx3DCrosstalkRL', p.trial.datapixx.crosstalk(:,1));    % ...only takes scalar gain param
    Datapixx('SetPropixx3DCrosstalkLR', p.trial.datapixx.crosstalk(:,end));    % ...only takes scalar gain param
    fprintLineBreak
    fprintf('Stereo Crosstalk correction implemented by Propixx firmware:\n');
    fprintf('\tL-(gain*R): [')
    fprintf('%05.2f, ', p.trial.datapixx.crosstalk(:,1).*100)
    fprintf('\b\b]%%\n')
    fprintf('\tR-(gain*L): [')
    fprintf('%05.2f, ', p.trial.datapixx.crosstalk(:,end).*100)
    fprintf('\b\b]%%\n')
end

%% Final status checks & updates before flip and return
% Retrieve/update extended Datapixx settings ("VideoStatus")
if p.trial.datapixx.use
    p.trial.datapixx.videoStatus = Datapixx('GetVideoStatus');
end

Screen('Flip', p.trial.display.ptr, 0);
