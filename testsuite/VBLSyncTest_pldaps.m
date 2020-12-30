function VBLSyncTest_pldaps(n, numifis, loadjitter, clearmode, stereo, flushpipe, synchronous, usedpixx, usevulkan, screenNumber)
% VBLSyncTest([n=600][, numifis=0][, loadjitter=0][, clearmode=0][, stereo=0][, flushpipe=0][, synchronous=0][, usedpixx=0][, usevulkan=0][, screenNumber=max])
%
%
% Tweaked from standard PTB code designed to utilize existing PLDAPS window (e.g. pre-initialized
% & configured exactly as in experiment)
% 
%
% Tests syncing of Psychtoolbox to the vertical retrace (VBL) and demonstrates
% how to implement the old Screen('WaitBlanking') behaviour with
% Screen('Flip')...
%
% This script provides a means to test, how well PTB synchronizes
% stimulus onset and execution of Matlab/Octave with the vertical retrace
% (also known as vertical blank or VBL) on your specific hardware setup.
%
% The script first opens a double-buffered fullscreen window. Then it
% performs a monitor calibration timing loop to estimate the real monitor refresh
% interval (aka IFI): While the formula ifi = 1.0 / Screen('NominalFramerate') will
% return an ifi that is close to the real IFI, it will be a little bit off
% from the real value. The reason is unavoidable jitter in the manufacturing of
% graphics cards internal clock circuits, as well as some drift and jitter
% caused by instabilities and change in the power supply and operating
% temperature of your machine. To be on the safe side, we use a timing-loop
% to compute the real IFI as an average of the IFI's of a number of consecutive
% monitor refresh intervals.
%
% After the calibration you'll see a simple animation: A flashing rectangle
% moving from the top-left corner of the screen to the bottom-right corner.
% During this animation, which lasts 'n' frames, it collects information
% about the timing behaviour.
%
% At the end, a couple of graphs are shown to you that should allow you to
% assess the timing behaviour of your setup.
%
% The following parameters can be changed in order to simulate different
% loads and stimulus presentation timings to assess PTB's behaviour under
% different conditions.
%
%
% PARAMETERS:
%
% n = Number of samples to take. E.g., n=1000 == Draw an animation
% consisting of 1000 frames.
%
%
%
% numifis = Number of monitor refresh intervals (IFIs) between flips:
% 0 == Flip at each vertical retrace: This is the old PTB 1.0.50 behaviour.
% Values of numifis>0 will cause Screen('Flip') to wait for 'numifis'
% monitor refresh intervals before flipping the back- and front buffers in
% sync with the vertical retrace.
%
% This would be roughly equivalent to the following snippet of code in the old
% MacOS9-PTB:
% ...
% Screen('WaitBlanking', windowPtr, numifis);
% Screen('CopyWindow', windowPtr, myOffscreenWindowwithStimulusPtr);
% ...
%
%
%
% loadjitter = Simulated load with a random duration between 0 ms
% and loadjitter monitor refresh intervals: We wait for the specified
% amount of time to simulate the execution of other Matlab-code, e.g.,
% KbChecks, GetMouse, Matlab calculations ...
%
%
%
% clearmode = Change the behaviour of Flip after flipping.
% clearmode = 0 will clear your stimulus drawing surface to background color after flip.
% This is the behaviour as found in PTB 1.0.50. After Flip you start with an
% empty image and can draw a completely new stim.
% 
% clearmode = 1 will not clear after a flip, but keep the contents of your stimulus
% image after the Flip: This allows you to incrementally update/draw stimuli.
%
% clearmode = 2 will neither clear nor keep the drawing surface after Flip, but leave the
% cleanup work to you. To be precise: The drawing surface will contain the
% stimulus image that was *just shown* on the screen. Think of Flip as if it would
% flip the front- and back-side of a sheet of paper: The current front side
% shows the stim to your subject, the current back side is where you draw.
% clearmode 2 will allow you to update the back side, which was the front
% side before the flip happened! This mode is useful if you want to save
% about 0.5-2 ms of time needed for mode 1 or 2 if you draw stimuli on very
% tight deadlines. 
%
%
% stereo = Test timing of display of stereoscopic stimuli.
% stereo = 0 will show you a standard monoscopic display.
% stereo = 1 will use the OS-X stereo output facilities to show stereoscopic
% stimuli: OS-X will quickly alternate between two images at each monitor refresh,
% one for the left-eye, one for the right-eye, while generating proper
% control signals for LCD shutter glasses. This should work with MacOS-X
% compatible stereo display hardware, e.g., CrystalEyes shutter glasses.
%
%
% flushpipe = Mark end of drawing commands to improve presentation timing.
% PTB knows a new command Screen('DrawingFinished') which, when properly used,
% will give PTB hints on how to optimize drawing of stimuli: This allows to draw
% more complex stimuli at higher monitor refesh intervals with reliable presentation
% timing. flushpipe enables/disables use of this new command.
%
% flushpipe = 0 Don't mark end of drawing commands.
% flushpipe = 1 Mark end of drawing commands to improve timing.
%
%
% synchronous = 0 Don't wait for drawing completion in Screen('DrawingFinished')
% synchronous = 1 Wait for completion - Useful for benchmarking and debugging,
% but degrades performance significantly in real experiments.
%
%
% usedpixx = 1 Use a DataPixx/ViewPixx/ProPixx device for external
% timestamping of stimulus onset, as a correctness test for Screen('Flip')
% timestamping. Disabled (0) by default.
% usedpixx = 2 Additionally correct for the clock skew between the computer
% and DataPixx device.
%
%
% 'usevulkan' If 1, try to use a Vulkan display backend instead of the
% OpenGL display backend. See 'help PsychVulkan'.
%
%
% screenNumber =  Use a screen other than the default (max) for testing .
%
%
% EXAMPLES:
%
% VBLSyncTest(1000, 0, 0.6, 0, 0, 1, 0) -- Render 1000 consecutive frames,
% flip at each retrace , pausing for 0.6 frame durations after each flip.
% Clear the framebuffer after flip, don't use stereo output, but use the new
% "DrawingFinished" command.
%
% VBLSyncTest(100, 10, 6, 1, 0, 1, 0) -- Render 100 frames,
% flip only every 10th monitor refresh interval ("WaitBlanking" for 10 refresh intervals),
% pause for 6 frame durations after each flip. Don't clear the framebuffer after
% flip, don't use stereo output, but use the new "DrawingFinished" command.
%
% VBLSyncTest(100, 10, 6, 1, 1, 1, 0) -- Render 100 frames,
% flip only every 10th monitor refresh interval ("WaitBlanking" for 10 refresh intervals),
% pause for 6 frame durations after each flip. Don't clear the framebuffer after
% flip, use stereo output via "frame-sequential stereo", use the new "DrawingFinished" command.
%
%
%
% THE PLOTS:
%
% Explanation of how to read the plots:
%
% Figure 1 shows the time delta between start of the VBL of successive
% Flip's: This value should be close to the requested delta as specified to
% Flip, e.g., you set numifis=0 or numifis=1 on a 100 Hz monitor. Then
% delta should be close to 1000 ms / 100 Hz = 10 ms. If numifis=2, it
% should be close to two monitor refresh intervals = 20 ms...
% 
% The green horizontal line denotes the proper delta value for your monitor
% refresh rate and 'numifis' value. The blue graph shows measured deltas. A
% jitter of less than +/- 1 ms indicates proper stimulus presentation timing -
% no skipped frames.
%
%
% Figure 2 shows the rasterbeam positions when flip took its internal
% timestamp: The values should be usually above the screen height, e.g., on
% a monitor resolution of 1200 x 1024 pixels, values should be above 1024.
% Values way below 1024 are also ok (e.g., < 100). This just means that
% your computer is either pretty slow, or connected to a flat-panel.
% Lots of randomly distributed values between 0 and 1024 would indicate sync trouble.
%
%
% Figure 3: Shows the estimated difference between requested presentation
% deadline and the real presentation deadline (start of VBL). Positive
% values indicate a deadline-miss and give you an indication of how much
% the deadline has been missed. Negative (or zero) values indicate that the
% deadline has been met. While the sign of this value is useful for assessing
% timing, the value itself is only meaningful for people who can read and
% fully understand the C source code and logic of 'Flips' implementation...
%
% Figure 4: Shows the difference (in milliseconds) between estimated
% start of VBL and return of the Flip command to Matlab. This is some
% indication of the processing overhead of OpenGL, the Operating system and
% Psychtoolbox when executing 'Flip'. It's also a lower bound for the
% timing delay when trying to synchronize start of acquisition devices to
% VBL.
%
% Figure 5: Shows the difference (in milliseconds) between estimated
% stimulus onset (aka end of vertical retrace, scanning beam starts
% at top of screen) and the end of Flip. This is the crucial value, if you
% want to sync something like sound-playback, triggering of some
% data-acquisition device (fMRI, MEG, EEG, ...) to stimulus onset. It
% should give you a feeling of how well you can sync. It's possible that
% negative values are reported on fast machines - Flip returns ahead of
% time (while monitor is still in retrace state). This is fine because your
% Matlab-Code for triggering something will add additional delays... Values
% should be below 1-2 milliseconds on reasonably modern and correctly
% configured hardware.
%
% Figure 6: Is only displayed when flag 'synchronous=1' This figure shows
% the total accumulated time for all drawing commands from the last Flip to the
% DrawingFinished command. It allows you to get a feeling on how hard the
% graphics hardware has to work for drawing your stim - and if it is
% possible at all to draw the stim on your hardware, given your time
% constraints. During a normal experiment (without sync-flag), the
% execution times of your Matlab code (as measured, e.g., by tic and toc)
% and of the drawing commands don't add up, because the graphics hardware
% works in parallel to the Matlab code.
%
%
% Please read the code of this M-File carefully as an example of how to get
% the best possible presentation timing on PTB-OSX.
%

%
%
% Date:   05/09/05
% Author: Mario Kleiner  (mario.kleiner.de@gmail.com)
%

%%% VBLSyncTest(1000, 0, 0.6, 0, 0, 1, 0)

if nargin < 1 || isempty(n)
    n = 600;
elseif isa(n, 'pldaps')
    p = n;
    w = p.trial.display.ptr;
    n = 5000;
    blockSize = 400;
    stereo = p.trial.display.stereoMode;
    white = p.trial.display.white;
    bgColor = p.trial.display.bgColor;
    ifi = p.trial.display.ifi;
    screenheight = diff(p.trial.display.winRect([2,4]));
    usedpixx = p.trial.datapixx.use * 2;
end

if ~exist('numifis','var') || isempty(numifis)
    numifis = 1;
end

if ~exist('loadjitter','var') || isempty(loadjitter)
    loadjitter = 0;
end

if ~exist('clearmode','var') || isempty(clearmode)
    clearmode = 0;
end

if ~exist('stereo','var') || isempty(stereo)
    stereo = 0;
end

if ~exist('flushpipe','var') || isempty(flushpipe)
    flushpipe = 0; % unused
end

if ~exist('synchronous','var') || isempty(synchronous)
    synchronous = 0; % unused
end

if ~exist('usedpixx','var') || isempty(usedpixx)
    usedpixx = 2;
end

if ~exist('usevulkan','var') || isempty(usevulkan)
    usevulkan = 0;
end

if ~exist('screenNumber','var')
    screenNumber = [];
end

if ~exist('blockSize','var') || isempty(blockSize)
    blockSize = 400;
end
if ~exist('white','var') || isempty(white)
    white = 2^8-1;
end
if ~exist('bgColor','var') || isempty(bgColor)
    bgColor = 127*[1 1 1];
end


    %PsychDefaultSetup(1);
%     RestrictKeysForKbCheck(KbName('ESCAPE'));
    
    if ~exist('w','var') || isempty(w)
        % Get the list of Screens and choose the one with the highest screen number.
        % Screen 0 is, by definition, the display with the menu bar. Often when
        % two monitors are connected the one without the menu bar is used as
        % the stimulus display.  Chosing the display with the highest display number is
        % a best guess about where you want the stimulus displayed.
        screens=Screen('Screens');
        if isempty(screenNumber)
            screenNumber=max(screens);
        end
        screensize=Screen('Rect', screenNumber);
        
        % Query size of screen:
        screenheight=screensize(4);
        
        % Open double-buffered window: Optionally enable stereo output if
        % stereo == 1.
        
        PsychImaging('PrepareConfiguration')
        if usedpixx
            % Use DataPixx for external timestamping for quick basic correctness
            % tests.
            PsychImaging('AddTask', 'General', 'UseDataPixx');
        end
        
        if usevulkan
            % Use PsychVulkan display backend instead of standard OpenGL:
            PsychImaging('AddTask', 'General', 'UseVulkanDisplay');
        end
        
        w=PsychImaging('OpenWindow',screenNumber, 0, [],[],[], stereo);
        
    end

    % Query effective stereo mode, as Screen() could have changed it behind our
    % back, e.g., if we asked for mode 1 but Screen() had to fallback to
    % mode 11:
    winfo = Screen('GetWindowInfo', w);
    stereo = winfo.StereoMode;
    
    if ~exist('ifi','var') || isempty(ifi)
        ifi = winfo.VideoRefreshFromBeamposition;
    end
    
    % Clear screen to black background color: If in stereo mode, we only
    % clear the left-eye buffer...
    Screen('SelectStereoDrawBuffer', w, 0);
    Screen('FillRect', w, bgColor);
    Screen('Flip', w);
    
    % Switch to realtime-priority to reduce timing jitter and interruptions
    % caused by other applications and the operating system itself:
    if IsLinux
        Priority(19);
    else
        Priority(MaxPriority(w));
    end
    
    % %     % Query nominal framerate as returned by Operating system:
    % %     % If OS returns 0, then we assume that we run on a flat-panel with
    % %     % fixed 60 Hz refresh interval.
    % %     framerate=Screen('NominalFramerate', w);
    % %     if (framerate==0)
    % %         framerate=60;
    % %     end
    % %
    % %     ifinominal=1 / framerate;
    % %     fprintf('The refresh interval reported by the operating system is %2.5f ms.\n', ifinominal*1000);

    % Start Screen timing for GPU rendering operations
    Screen('GetWindowInfo', w, 5);

    % Perform a calibration loop to determine the "real" interframe interval
    % for the given gfx-card + monitor combination:
    Screen('TextSize', w, 24);
    Screen('DrawText', w, 'Measuring monitor refresh interval... This can take up to 20 seconds...', 10, 10, 255);
    
    if (stereo>0)
        % Show something for the right eye as well in stereo mode:
        Screen('SelectStereoDrawBuffer', w, 1);
        Screen('FillRect', w, bgColor);
        Screen('DrawText', w, 'Stereo yeah!!!', 10, 40, 255);
    end
    
    
    % %     % Measure monitor refresh interval again, just for fun...
    % %     % This will trigger a calibration loop of minimum 100 valid samples and return the
    % %     % estimated ifi in 'ifi': We require an accuracy of 0.1 ms == 0.0001
    % %     % secs. If this level of accuracy can't be reached, we time out after
    % %     % 20 seconds...
    % %     %[ ifi nvalid stddev ]= Screen('GetFlipInterval', w, 100, 0.0001, 5);
    % %     [ ifi, nvalid, stddev ]= Screen('GetFlipInterval', w);
    % %     fprintf('Measured refresh interval, as reported by "GetFlipInterval" is %2.5f ms. (nsamples = %i, stddev = %2.5f ms)\n', ifi*1000, nvalid, stddev*1000);
    
    % Init data-collection arrays for collection of n samples:
    ts=zeros(1,n);
    beampos=ts;
    missest=ts;
    flipfin=ts;
    dpixxdelay=ts;
    td=ts;
    so=ts;
    sodpixx = ts;
    boxTime = ts;
    
    frameRenderTime = ts;    
    
    % Compute random load distribution for provided loadjitter value:
    if isscalar(loadjitter)
        % as multiple of frame ifi
        wt=rand(1,n)*(loadjitter*ifi);
    else
        % as limits (in ms) of a uniform distribution
        wt = (rand(1,n)*diff(loadjitter)+loadjitter(1)) / 1000;
    end


    % Perform some initial Flip to get us in sync with retrace:
    % tvbl is the timestamp (system time in seconds) when the retrace
    % started. We need it as a reference value for our WaitBlanking
    % emulation:
    tvbl=Screen('Flip', w);
    
    % Test-loop: Collects n samples.
    for i=1:n
        % Presentation time calculation for waiting 'numifis' monitor refresh
        % intervals before flipping front- and backbuffer:
        % This formula emulates the old PTB-MacOS9 Screen('WaitBlanking', numifis)
        % behaviour as closely as possible.
        % The 'Flip' command takes a presentation timestamp 'tdeadline' as
        % optional argument: If tdeadline == 0 or is left out, Flip will
        % flip at the next possible retrace (PTB 1.0.5 and earlier behaviour).
        % If tdeadline is > 0, then Flip will wait until the system time
        % 'tdeadline' is reached and then flip the buffers
        % at the next possible VBL. This allows to specify absolute points in time
        % at which flip should occur. If you want the old behaviour of
        % Screen('WaitBlanking', w, numifis) back, then just calculate a
        % proper presentation timestamp "tdeadline" relative to the time of last flip
        % 'tvbl', as demonstrated here:        
        tdeadline=tvbl + numifis * ifi - 0.5 * ifi;
        
        % If numifis == 0, flip on next retrace. This should be the same
        % as numifis == 1, but might make a difference in robustness if
        % the stimulus is **very** complex and the load for the system is at
        % the limit that it can handle in a single video-refresh interval.
        if numifis==0
            % If user supplied numifis=0, we force tdeadline=0, so Flip
            % will actually ignore the deadline and just Flip at the next
            % possible retrace...
            tdeadline=0;
        end
        
        if usedpixx
            % Ask for a Datapixx onset timestamp for next 'Flip':
            PsychDataPixx('LogOnsetTimestamps', 1);
        end
        
        % Flip: The clearmode argument specifies if flip should clear the
        % drawing buffer after flip (=0 - default), keep it "as is"
        % for incremental drawing/updating of stims (=1) or don't do
        % anything to the framebuffer at all (=2).
        % We return the timestamp, when VBL starts in tvbl: This is when
        % the front- and back drawing surfaces get exchanged and it is the
        % crucial reference value for computing the 'tdeadline'
        % presentation deadline for the next 'Flip' command.
        % The rasterbeam-position (scanline) when the measurement was taken is returned in beampos(i),
        % the time when flip returned to Matlab is returned in flipfin(i),
        % estimated stimulus onset time aka end of VBL is returned in so(i).
        %
        % The first value "tvbl" is needed for tdeadline calculation if
        % one wants to emulate WaitBlanking of old PTB - see formula above.
        % beampos > screen height means that flip returned during the VBL
        % interval. Small values << screen height are also ok,
        % they just indicate either a slower machine or some types of flat-panels...

        [ tvbl, so(i), flipfin(i), missest(i), beampos(i)]=Screen('Flip', w, tdeadline, clearmode);

        if usedpixx
            % Ask for a Datapixx onset timestamp from last 'Flip':
            [boxTime(i), sodpixx(i)] = PsychDataPixx('GetLastOnsetTimestamp');
            dpixxdelay(i) = GetSecs;
        end
        
        winfo = Screen('GetWindowInfo', w);
        frameRenderTime(i) = winfo.GPULastFrameRenderTime;
        
        % Start Screen timing for GPU rendering operations
        if i<n
            Screen('GetWindowInfo', w, 5);
        end
        % Record timestamp for later use:
        ts(i) = tvbl;
        
        % Draw some simple stim for next frame of animation: We draw a
        % simple flashing rectangle that moves over the screen. The same
        % rectangle is drawn with some offset for the right-eye if stereo
        % display is requested:
        Screen('SelectStereoDrawBuffer', w, 0);
        pos=mod(i, screenheight-blockSize);
        blockPos = pos + blockSize*[.05, .05, 1, 1];
        blockCol = mod(GetSecs, white);
        Screen('FillRect', w, blockCol, blockPos);
        %         Screen('FillRect', w, mod(i, 255), [pos+20 pos+20 pos+400 pos+400]);
        if (stereo>0)
            % Show something for the right eye as well in stereo mode:
            Screen('SelectStereoDrawBuffer', w, 1);
            Screen('FillRect', w, blockCol, blockPos + blockSize*[.05 0 .05 0])
            %             Screen('FillRect', w, mod(i, 255), [pos+40 pos+20 pos+420 pos+400]);
        end

        % ALWAYS use 'DrawingFinished'
        Screen('DrawingFinished', w, clearmode);
        % shift all render timing estimates to proper Screen('GetWindowInfo'...) method
        
        % Sleep a random amount of time, just to simulate some work being
        % done in the Matlab loop:
        WaitSecs(wt(i));
        % And give user a chance to abort the test by pressing any key:
        if KbCheck
            break;
        end
    end % Draw next frame...

    % calculate clock skew corrected Datapixx onset timestamps
    if usedpixx>1
        sodpixx = PsychDataPixx('BoxsecsToGetsecs', boxTime);
    end

    RestrictKeysForKbCheck([]);
    % final Flip
    Screen('Flip', w, 0, clearmode);
    
    % Shutdown realtime scheduling:
%     Priority(0)
% 
%     % Close display: If we skipped/missed any presentation deadline during
%     % Flip, Psychtoolbox will automatically display some warning message on the Matlab
%     % console:
%     sca;

    % Restrict to actual number of collected samples:
    % - skip first sample
    n = i;
    ii = 2:n;
    ts = [nan, ts(ii)];
    so = [nan, so(ii)];
    flipfin = [nan, flipfin(ii)];
    missest = [nan, missest(ii)];
    beampos = [nan, beampos(ii)];
    td = [nan, td(ii)];
    dpixxdelay = [nan, dpixxdelay(ii)];
    sodpixx = [nan, sodpixx(ii)];
    boxTime = [nan, boxTime(ii)];

    frameRenderTime = [nan, frameRenderTime(ii)];

    % Count and output number of missed flip on VBL deadlines:
    numbermisses=0;
    numberearly=0;

    if numifis > 0
        if (stereo == 11) && (numifis == 1)
            % Special case: Stereomode 11 can't do better than one swap
            % every two refresh cycles, so take this into account:
            ifi = ifi * 2;
        end

        for i=2:n
            if (ts(i)-ts(i-1) > ifi*(numifis+0.5))
                numbermisses=numbermisses+1;
            end

            if (ts(i)-ts(i-1) < ifi*(numifis-0.5))
                numberearly=numberearly+1;
            end
        end
    else
        if stereo == 11
            % Special case: Stereomode 11 can't do better than one swap
            % every two refresh cycles at best, so take this into account:
            ifi = ifi * 2;
        end

        for i=2:n
            if (ts(i)-ts(i-1) > ifi*1.5)
                numbermisses=numbermisses+1;
            end

            if (ts(i)-ts(i-1) < ifi*(numifis-0.5))
                numberearly=numberearly+1;
            end
        end
    end

    % Plot all our measurement results:

    % Figure 1 shows time deltas between successive flips in milliseconds:
    % This should equal the product numifis * ifi:
    figN = 1;
    figure
    hold on
    plot((ts(2:n) - ts(1:n-1)) * 1000);
    ni = numifis;
    if numifis < 1
        ni = 1;
    end
    
    if (ni < 2 ) && (stereo == 11)
        % Special case: Stereomode 11 can't do better than one swap
        % every two refresh cycles, so take this into account:
        ni = 2;
    end
    
    plot(ones(1,n)*ifi*ni*1000);
    tiStr = ('Delta between successive Flips in milliseconds:');
    title({sprintf('Figure %d',figN), tiStr});
    box off
    hold off
    
    % Figure 2 shows the recorded beam positions:
    if winfo.VBLEndline > -1
        hasbeampos = 1;
    else
        hasbeampos = 0;
    end

    if hasbeampos
        figure
        plot(beampos);
        tiStr = ('Rasterbeam position when timestamp was taken (in scanlines):');
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off

    end

    if (numbermisses > 1) || (numberearly > 0)
        % Figure 3 shows estimated size of presentation deadline-miss in
        % milliseconds:
        figure
        hold on
        plot(missest*1000);
        plot(zeros(1,n));
        tiStr = ('Estimate of missed deadlines in milliseconds (negative == no miss):');
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off
        hold off
    end

    if isequal(ts, so)
        % Same info in vbltime and stimulus onset time. Only
        % do one plot and label it in a less confusing manner:
        figure
        plot((flipfin - so)*1000);

        if IsLinux && (Screen('Preference', 'VBLTimestampingmode') == 4)
            % Linux mode 4: OpenML or Wayland presentation feedback, so is stimulus onset:
            tiStr = ('Time delta between stimulus onset and return of Flip in milliseconds:');
        elseif (IsLinux || IsOSX) && ismember(Screen('Preference', 'VBLTimestampingmode'), [1, 3]) && (winfo.VBLCount > 0)
            % Linux or OSX, vbl timestamping requested and apparently working. so is vblank time:
            tiStr = ('Time delta between start of VBL and return of Flip in milliseconds:');
        else
            % Windows or other os'es without working high precision timestamping: so is raw timestamp:
            tiStr = ('Time delta between return from swap completion and return of Flip in milliseconds:');
        end
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off
        
    else
        % Figure 4 shows difference in ms between finish of Flip and estimated
        % start of VBL time:
        figure
        plot((flipfin - ts)*1000);
        tiStr = ('Time delta between start of VBL and return of Flip in milliseconds:');
        figN = figN+1;
        title({sprintf('Figure %d-a',figN), tiStr});
        box off

        % Figure 5 shows difference in ms between finish of Flip and estimated
        % stimulus-onset:
        figure
        plot((flipfin - so)*1000);
        tiStr = ('Time delta between stimulus onset and return of Flip in milliseconds:');
        title({sprintf('Figure %d-b',figN), tiStr});
        box off
    end

    % Figure 6 shows duration of drawing commands (rendering times)
    % when calling "DrawingFinished"
        figure
        plot(frameRenderTime * 1000);
        tiStr = ('Frame render duration in milliseconds:');
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off
        
    if usedpixx
        figure;
        plot((so - sodpixx) * 1000);
        tiStr = ('Time delta in msecs onset according to Flip - onset according to DataPixx:');
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off
        fprintf('Average discrepancy between Flip timestamping and DataPixx is %f msecs, stddev = %f msecs.\n', mean((so - sodpixx) * 1000), std((so - sodpixx) * 1000));

        figure;
        plot((dpixxdelay - so) * 1000);
        tiStr = ('Time delta between stimulus onset and return of Datapixx timestamping in milliseconds:');
        figN = figN+1;
        title({sprintf('Figure %d',figN), tiStr});
        box off
    end

    % Output some summary and say goodbye...
    fprintf('PTB missed %i out of %i stimulus presentation deadlines.\n', numbermisses, n);
    fprintf('One missed deadline is ok and an artifact of the measurement.\n');
    
    fprintf('PTB completed %i stimulus presentations before the requested target time.\n', numberearly);
    if numberearly > 0
        fprintf('CAUTION: Completing flips too early should *never ever happen*! Your system has\n');
        fprintf('CAUTION: a serious bug or misconfiguration in its graphics driver!!!\n');
    end
    fprintf('Have a look at the plots for more details...\n');
    
    % Done.

return
