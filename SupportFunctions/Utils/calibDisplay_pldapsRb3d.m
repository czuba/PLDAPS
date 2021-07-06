function [simpleGamma, measInt, measCie, lvlSet] = calibDisplay_pldapsRb3d(p, mode)
% function [simpleGamma, measInt, measCie, lvlSet] = calibDisplay_pldapsRb3d(p, mode)
%  [p] = currently running PLDAPS structure
%  [mode]: 0 = grey[k] only, no SPD;
%          [1 2 3 4] == [R G B K], luminance, xyCIE, & SPD
% 
% Controls PR655 to measure luminance at range of intensity values
% from one or multiple locations on the screen. 
% Should work on all PLDAPS setups, regardless of overlay settings
% 
%   ** uses tweaked subfunction versions of PR655init.m and parseSpd.m 
% 
% USAGE: 
%  Must be run after pausing a running PLDAPS session.
%  - This is a feature, not a bug, because it ensures that measurements are made under
%    the exact same setup conditions as your experiment.
%  This is a special version for use with ProPixx RB3D stereo projector modes
%  - includes some extra setup & measurements for monocular calibration,
%    & binocular crosstalk correction
%  - ...still needs [more] commenting & explanation on how to properly use/implement (TBC 2021)
% 
%  When making initial calibration measurements, its important to be sure that 
%  all other/prior calibrations are disabled when the new measurements are made.
%  - e.g. linearization via gamma table or power, or binocular crosstalk corrections
%  - Otherwise the result will be confounded by [unknown] prior display corrections
%  - Easiest way to do this is to make sure the .display.gamma.power field is empty
%    before your code executes p.run
%  - In the modularDemo, you would want to set:  pss.display.gamma.power = [];
%  
%  Save the returned greyscale [simpleGamma] value in your rigPrefs as:
%       .display.gamma.power
% 
%  AFTER you've run your calibration & added the results to your rigPrefs,
%  run this calibration function again to validate. This will help ensure that your
%  validation is performed in EXACTLY the same state as your experiment.
% 
% 
% EXAMPLE:
%   p = yourPldapsExperiment;
%   % [ press 'p' key to pause PLDAPS ]
% 
%   simpleGamma = calibDisplay_pldaps(p);
%   % [ follow command window prompts for calibration ]
%   % [ calibration will return to pause state when finished ]
% 
%   % copy return argument to Matlab 'base' workspace
%   assignin('base', 'simpleGamma', simpleGamma)
% 
%   % [ resume & exit experiment gracefully (Type "dbcont" to continue..., then press 'q' to quit) ]
% 
%   % from the command window (after completely exiting the PLDAPS experiment),
%   % add [simpleGamma] to your PLDAPS rigPrefs
%   myGamma.display.gamma.power = simpleGamma;
%   createRigPrefs(myGamma)
%   % [ follow command window instructions to apply & save new rigPrefs values ]
%   
%   % Done! Future experiments on this rig will apply this gamma correction automatically
%   % Rerun the calibration measurements with the new settings applied to validate your results
% ---
% 
% 
% ------------------
% TBC 08-05-2008   Thad Czuba wrote it (UT Austin)
% TBC 09-13-2016   Scruffy update. Functionified & subfunctioned dependencies.
% 2017-10-18  TBC  PLDAPS update
% 2019-02-12  TBC  Updated version for Propixx RB3D & crosstalk calibration measurements
% 2021-03-12  TBC  General refinements & cleanup
%                  Updated default save dir to ./PLDAPScalibrate on same level as PLDAPS installation
%                  e.g.  ~/MLtoolbox/PLDAPScalibrate
% 


%% Parse inputs & check compatibility
% !!ONLY!! compatible with Rb3d stereo calibration
if ~isfield(p.trial.datapixx,'rb3d') || ~p.trial.datapixx.rb3d
    error('RB3d mode needs special treatment to control LED illuminator\nduring R/G/B calibration & crosstalk measurements.\nYou appear to have RB3D disabled, so this is probably not the function for  you...\n')
end

if nargin<2 || isempty(mode)
    mode = 0;
end


%% Paths & Filenames
cal.dispname = input('Input name of display (i.e. rig name): ','s');
cal.filename = [cal.dispname, datestr(p.trial.session.initTime, 'yyyymmdd')];

% Default save dir:
% - Same enclosing* directory as your PLDAPS installation (presumably git repo; likely  ~/MLtoolbox )
%  - **but not inside PLDAPS**, else could accidentally delete all calibration data when fixing/replacing your repo
%  - <PLDAPS root dir>/../PLDAPScalibrate/<YYYY>
cal.calibDir = fullfile( fileparts(p.trial.pldaps.dirs.proot), 'PLDAPScalibrate', datestr(date,'yyyy'));
if ~exist(cal.calibDir,'dir')
    mkdir(cal.calibDir);
end

cal.rawDir = fullfile(cal.calibDir, 'Raw');
if ~exist(cal.rawDir,'dir')
    mkdir(cal.rawDir);
end
cal.src         = fullfile(cal.calibDir,    [cal.filename,'.mat']);
cal.srcWkspc    = fullfile(cal.rawDir,      [cal.filename,'wkspc.mat']);


%% Defaults
leaveRoomTime = 5; % sec
cal.measmin = BlackIndex(p.trial.display.ptr); % 0;  % min intensity value
cal.measmax = WhiteIndex(p.trial.display.ptr); % 1;  % max intensity value
defaultScreen = max(Screen('Screens'));
defaultMeasRes = 2^8; % number intensity levels
defaultSampDur = 1200; % sample duration in ms
defaultAvg = 2; % number of repeats through intensity range to average: 1 for speed, 2-3 recommended
% ...distinct from the averaging done within photometer because this timescale will address
% slow power fluctuations that can occur over course of a long calibration (e.g. 10+bit measRes)

% sample grid (Locations of measurement points & selected points to measure)
defaultGridSize = [3, 1]; % n-grid locations in [x y] screen dimensions
defaultGridPts = [2, 2];%[2];%[1:3];

% sync photometer to display refresh rate
% - this can be buggy/problematic, esp for slow LCD displays or exotic refresh paradigms
% - ALWAYS disabled for RB3D (doesn't work with 'galloping' timing imposed by Propixx RB3D mode)
defaultSyncToSource = 1;

dpxMaskIndex = [6, 5, 3, 0];   %[6, 5, 3, 0] == [R G B K];    %[0] == [K];   %
if ~mode
    % greyscale only
    dpxMaskIndex = dpxMaskIndex(4);
else
    dpxMaskIndex = dpxMaskIndex(mode);
end

% List defaults
fprintf('Default settings:\n\tScreen:\t\t%g\n\tGrid Size:\t%gx%g\n\tGrid Locs:\t%s\n\tInt [Min,Max]:\t[%g,%g]\n\tMeas Res:\t%g\n\tSample Avg:\t%g\n',...
    defaultScreen,defaultGridSize,num2str(defaultGridPts),[cal.measmin,cal.measmax],defaultMeasRes,defaultAvg);   %,defaultDist)

goDef = input('Accept all defaults? [1=yes, 0=no]');
if ~goDef
    % Grid size
    gridSize = input(sprintf('Size of grid dimensions (format: [x,y]; def: [%d,%d]): ', defaultGridSize));
    if isscalar(gridSize),      gridSize = gridSize*[1 1];
    elseif isempty(gridSize),   gridSize = defaultGridSize; end
    disp(eval(mat2str(reshape(1:prod(gridSize),gridSize))))
    % ensure defaultGridPts is within current range of gridSize
    if any(defaultGridPts > prod(gridSize)), defaultGridPts = round(prod(gridSize)/2);   end
    
    % Grid points to use
    gridPts = input(sprintf('Grid locations to sample (available: %s; def: %s): ', mat2str(1:prod(gridSize)), mat2str(defaultGridPts)));
    
    %  What intensity levels do you want to measure?
    measRes = input(sprintf('How many intensities do you want to measure? (min:16, def:%g) ',defaultMeasRes));
    
    % Sync to source refresh rate?
    syncToSource = input(sprintf('Synchronize measurements to source refresh rate? [1=yes, 0=no; def: %d]: ', defaultSyncToSource));
    
    % Duration of each sample?
    sampDur = input(sprintf('Duration of each sample exposure [%d ms, int]: ', defaultSampDur));
    
    % Samples to average
    nAverage = input(sprintf('Num of samples to average [%d]: ', defaultAvg));
    
end

% initGamma = repmat(linspace(0,1,256)',1,3);
if ~exist('gridSize','var') || isempty(gridSize),    gridSize = defaultGridSize; end

if ~exist('gridPts','var') || isempty(gridPts),    gridPts = defaultGridPts; end

if ~exist('measRes','var') || isempty(measRes),    measRes = defaultMeasRes; end
% elseif measRes>2^ScreenDacBits(whichScreen), disp('Cannot exceed bit resolution. Exiting...'),return,    end

if ~exist('syncToSource','var') || isempty(syncToSource),   syncToSource = 1;  end
if isfield(p.trial.display,'rb3d') && p.trial.display.rb3d
    syncToSource = 0; % ALWAYS disabled for RB3D (doesn't work with 'galloping' timing imposed by Propixx RB3D mode)
end

if ~exist('sampDur','var') || isempty(sampDur),    sampDur = defaultSampDur; end

if ~exist('nAverage','var') || isempty(nAverage),   nAverage = defaultAvg;  end


%% Establish variables
cal.measRes = measRes;
cal.syncToSource = syncToSource;

% Update some general system & setup variables
winPtr = p.trial.display.ptr;
winRect = p.trial.display.winRect;
% % ensure desired led calibration setting, and update datapixx video status
% Datapixx('setgrayledcurrents',0);   pause(.1)
cal.winInfo = Screen('GetWindowInfo', p.trial.display.ptr);
if p.trial.datapixx.use
    % ensure desired led calibration setting, and update datapixx video status
    ledCurrentSetIndex = 0;
    Datapixx('setgrayledcurrents', ledCurrentSetIndex);
    if ~ledCurrentSetIndex
        fprintf('~~~\tDefault Propixx LED current set selected:\tledCurrentSetIndex = %s\n~~~\tUse this mode for ACQUIRING baseline calibration measurements.\n', ledCurrentSetIndex)
    else
        fprintf(2, '~~~\tCalibrated Propixx LED current set selected:\tledCurrentSetIndex = %s\n~~~\tUse this mode for TESTING calibration accuracy.\n~~~\t(LED currents must be set manually in vputil app. ...see calibration readme or vpixx documentation for instructions.)\n')
    end
    pause(.4)
    cal.dpxInfo = Datapixx('GetVideoStatus');
end
cal.gridSize = gridSize;
cal.gridPts = gridPts;

% Set up sample area grid.
% (krufty, but didn't want to recode...)
xmax = gridSize(1);
ymax = gridSize(2);
[yy, xx] = meshgrid(1:ymax, 1:xmax);

grid = cell(gridSize); %z = 1;
for i = 1:numel(grid) %prod(gridSize);
    grid{i} = [xx(i)/xmax, yy(i)/ymax];
end

% Make grid of test circle coords
stimsz = min(grid{1}.*winRect(3:4)-10);
stimdiff = (grid{1}.*winRect(3:4)-stimsz)/2;
circRect = [-(stimsz+stimdiff(1)+5) -(stimsz+stimdiff(2)+5) -stimdiff(1) -stimdiff(2)];

% Show all spots
Screen('FillRect', winPtr, [0 0 0]')
for i = 1:numel(grid)
    circ = repmat(grid{i}.*winRect(3:4),1,2)+circRect;
    Screen('FillOval', winPtr, cal.measmax*0.7*[1 1 1]', circ)
end
t=Screen('Flip', winPtr);



%% Initialize PR655 in Remote Mode
disp('Initializing PR-655')
rm = PR655init;
disp(rm);

pdur = 0.3; % wait for photometer to process input (sec)
pause(pdur)

% turn photometer backlight off
PR655write('B00');  % this setting doesn't always 'stick', suggest draping with cloth during recordings, just to be safe
pause(pdur)

% duration of measurements in msec
% PR655write('SE0');  disp('Adaptive exposure timing...');
% % NOTE: Adaptive duration not recommended for luminance calibration; it will fail when very dim at low end
PR655write(sprintf('SE%d', sampDur));    % PR655write('SE2400'); % PR655write('SE600');
pause(pdur)

% number of measurements to average
PR655write('SN3');  % average 3 samples within photometer
pause(pdur)

% Sync to source frequency
% - this can be buggy/problematic, esp for slow LCD displays or exotic refresh paradigms
% - ALWAYS disabled for RB3D (doesn't work with 'galloping' timing imposed by Propixx RB3D mode)
if cal.syncToSource
    syncFreq = PR655getsyncfreq;
    if ~isempty(syncFreq) && syncFreq ~= 0
        fprintf('Measured sync frequency of: %3.0fHz.\nSyncing PR655 to source w/ adaptive exposure timing...\n',syncFreq);
        PR655write('SE0');
        pause(pdur)
        PR655write('SS1');
    else
        cal.syncToSource = 0;
        PR655write('SS0');
        disp('Warning: Could not sync photometer to source frequency.');
    end
else
    PR655write('SS0');
end
pause(pdur)
pause(pdur)

% Read out & confirm photometer settings
PR655write('D602');
cal.pmeterSettings = PR655read;
disp(cal.pmeterSettings);

hasBegun = 0;
pause(2)

%% Begin calibration
% Baseline set of requested intensities
lvlSet = linspace(cal.measmin, cal.measmax, measRes);
cal.lvlSet = lvlSet;
% RB3D only works in grayscale mode, so we need to use LED mask to measure RGB outputs individually
% NOTE:  because of this, we also don't want/need to constrain intensity levels to a single channel

% Mask for measuring grey, R, G, B luminance outputs
cal.colmask = dpxMaskIndex;

for maski = circshift(1:length(dpxMaskIndex), 1) % if multiple, do grey first (but don't mess w/indexing)
    Datapixx('SetPropixxLedMask', dpxMaskIndex(maski));
    Datapixx('RegWrRd');
    
    for r = 1:length(gridPts)   % for each calibration location
        % Co-opt normal mechanism for measuring multiple locations on screen, to different eye stimulation conditions
        if r==1
            % left eye
            iBuffer = 0;
        elseif r==2
            % right eye
            iBuffer = 1;
        elseif r==3
            % both eyes
            iBuffer = [0,1];
        else
            % not defined for >3
            break
        end
        
        circ = repmat(grid{gridPts(r)}.*winRect(3:4),1,2)+circRect;
        
        % Pause for setup if first pass, or measuring multiple/distinct screen locations
        if ~hasBegun || length(unique(gridPts))>1
            Screen('FillRect',winPtr,p.trial.display.bgColor)
            Screen('FillOval',winPtr,cal.measmax*[1 1 1],circ)
            Screen('DrawLine',winPtr, 0.1*cal.measmax*[1 1 1],circ(1),circ(2),circ(3),circ(4),2);
            Screen('DrawLine',winPtr, 0.1*cal.measmax*[1 1 1],circ(1),circ(4),circ(3),circ(2),2);
            t = Screen('Flip',winPtr);
            
            disp('Aim Photometer at center of spot and press space bar.')
            fprintf('You will have %2.0f sec to leave.', leaveRoomTime)
            while 1
                [~, ~, key] = KbCheck(-1);
                if key(p.trial.keyboard.codes.spaceKey), break, end
            end
            hasBegun = 1;
            Beeper, Beeper
            pause(leaveRoomTime)
        end
        
        tic
        
        for a = 1:nAverage % repeated measurements
            % randomize order
            jnk =  randperm(length(lvlSet)); %1:length(lvl); fprintf('\n\nDEBUG: not randomizing order!\n\n'); %

            for i = 1:length(jnk) % randomized walk through intensity lvls
                lvlIdx = jnk(i);
                % Intensity level for this pass
                lvl = lvlSet(lvlIdx) * [1 1 1]';
                
                % Fill background when appropriate (but ALWAYS in both eyes)
                for buf = p.trial.display.bufferIdx
                    Screen('SelectStereoDrawBuffer', winPtr, buf);
                    if dpxMaskIndex(maski)~=0
                        % no grey surround (LEDmask makes colored surround, skewing measurement)
                        Screen('FillRect',winPtr, [0 0 0]);
                    else
                        % grey background
                        Screen('FillRect', winPtr, p.trial.display.bgColor);
                        % black target circle
                        Screen('FillOval', winPtr, [0 0 0], circ);
                    end
                end
                
                % Show circle of given intensity in appropriate eye(s)
                for buf = iBuffer
                    Screen('SelectStereoDrawBuffer', winPtr, buf);
                    Screen('FillOval', winPtr, lvl, circ);
                end
                
                % Flip screen
                t = Screen('Flip', winPtr);
                
                % Measure intensity
                PR655write('M1');
                StartTime = GetSecs;
                tmp = [];
                while isempty(tmp)
                    tmp = PR655read;
                    
                    if GetSecs-StartTime > 60
                        % if no measurement returned w/in 60 sec, save data collected & exit calibration
                        measError_saveAndReturn();
                    end
                end
                % retval is raw 1-by-5 values returned by photometer: [status, units, photometric brightness, CIE 1931 x, CIE 1931 y] (...per PR655 manual)
                retval{r, maski}(lvlIdx,:,a) = str2num(tmp);
                
                if a==1 && any(mode)
                    % Lets measure SPD while we're here too
                    PR655write('M5');
                    StartTime = GetSecs;
                    tmp = [];
                    while isempty(tmp)
                        tmp = PR655read;
                        
                        if GetSecs-StartTime > 60
                        % if no measurement returned w/in 60 sec, save data collected & exit calibration
                        measError_saveAndReturn();
                        end
                    end
                    % parseSpd == subfunction based on PR655parsespdstr (but w/o questionable conversion&splining)
                    [measSpd{r, maski}(:,lvlIdx), spdWfs] = parseSpd(tmp);           %#ok<*NASGU,*AGROW,*ASGLU>
                end
                % show vals in command window
                fprintf('\n');    fprintf([repmat('\b',[1,60]),'\n']); % clear previous printout
                fprintf('{%d, %d} lvl(%03.0f): [%05.1f, x %03.3f, y %03.3f]  %04.1f%% complete', r,maski,lvlIdx, retval{r, maski}(lvlIdx, 3:5, a), (i+(a-1)*measRes)/(measRes*nAverage)*100);
                
            end
            
        end
        
        % Parse & average repetitions from raw output
        %       [status, units, photometric brightness, CIE 1931 x, CIE 1931 y]
        measInt = cellfun(@(x) mean(x(:,3,:), 3), retval, 'uni',0);     % luminance
        measCie = cellfun(@(x) mean(x(:,4:5,:), 3), retval, 'uni',0);   % chromaticity
        
        % save alot (jic)
        try
            save(cal.srcWkspc);
        end
        
    end
end

if ~exist('measSpd','var')
    measSpd = [];
end


%% Fit the data
[fx, fo] = deal(nan(size(measInt)));
for i = 1:numel(measInt)
    measNorm = measInt{i}./max(measInt{i});
    [~, tmp] = FitGamma(lvlSet(:), measNorm, lvlSet(:), 2);
    fx(i) = tmp(1);
    fo(i) = tmp(2);
end
simpleGamma = 1./fx;

gam.src = cal.src;
gam.power = simpleGamma(end); % grey should always be last
gam.simpleGamma = simpleGamma;
gam.offset = fo;

% Save variables of interest and update full workspace file
try
    save(cal.src, 'measInt','measCie','measSpd','gam','cal');
    save(cal.srcWkspc);
    fprintf(2, '\tGamma correction & calib data saved to\n\t\t%s\n\n', cal.src);
    disp(eval(mat2str(simpleGamma)));
catch
    warning('Error occurred while trying to save calibration outputs.\n\tStopping w/in calibration fxn so you have a chance to manually save/inspect.\n',[]);   %#ok<CTPCT>
    keyboard
end




%% Clean up and save
PR655close;
%%

Screen('FillRect', winPtr, p.trial.display.bgColor)
for i = 1:numel(grid)
    circ = repmat(grid{i}.*winRect(3:4),1,2)+circRect;
    Screen('FillOval',winPtr,.8*[1,1,1],circ)
end
Screen('Flip',winPtr);


disp('Done.')

% % if ~isempty(emailfinish)
% %     try
% % 	send_mail(emailfinish,'Calibration complete',sprintf(['Congratulations, The Flying Spaghetti Monster has granted your wish and completed your calibration of ',...
% %         cal.dispname,'. \n\nThe data file and related gamma/conversion files are located at:\n\t',gammaDir,gamma.fname,'\n\n\nGoodbye,\n\tFSM']));
% %     catch, disp('done.')
% %     end
% % end


% % % % % % % % % % % % % % % %
%% Nested functions
%  (shares workspace with main function)
% % % % % % % % % % % % % % % %
    function measError_saveAndReturn
        % Breakout function for saving error data when measurement problems arise (no meas returned for 60 sec)
        fprintf(2, '\n\nNo measurement was returned...something broke. \t[%s]\n', datestr(now));
        errFile = fullfile(cal.rawDir, [cal.filename,'_errorWkspc_',datestr(now,'HHMMSS'),'.mat']);
        fprintf(2, 'Saving data to: \n\t%s\n', errFile);
        save( errFile );
        % Attempt to close photometer connection
        try
            PR655close;
        end
        return
        
    end % measError_saveAndReturn


end %main function

% % % % % % % % % % % % % % % %
% End of primary code block
% % % % % % % % % % % % % % % %




% % % % % % % % % % % % % % % %
%% Subfunctions
%  simplify dependencies (...lots of copy pasta here)
% % % % % % % % % % % % % % % %

function [spd, wfs] = parseSpd(readStr)
% spd = PR655parsespdstr(readStr,S)
%
% Parse the spectral power distribution string
% returned by the PR655. ...start at 380nm (standard), and parse vals until fail.
%
% [spd, wfs] == [rawPower, wavelength]
%
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655,
%                   removing any conversions, returns exactly what the PR puts out.
%

k = 0;
start = strfind(readStr,'380,');
while 1
    try
        k = k+1;
        wfs(k) = str2num(readStr(start+16*(k-1):start+2+16*(k-1)));
        spd(k) = str2num(readStr(start+4+16*(k-1):start+4+9+16*(k-1)));
    catch
        break
    end
end

end
% % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % %


%% PR655init
function retval = PR655init
%
% function retval = PR655init
%
% Simpler, more robust version of PR655init.m  Brief pauses after IO seems to improve reliability
% with fast/modern matlab communication (..>=2016).
%   ** This might not work on Windows machines. (...that much we have in common)
%   --TBC (circa 2017)
% 
% If still having trouble connecting, try a different USB port/board. Seems 3rd party usb expansion
% boards can mess with USB connectivity through IOPort (after unknown trigger conditions).
% 
% % Initialize serial port for talking to colorimeter.
% Returns whatever character is sent by colorimeter
%
% Per PhotoResearch, handshaking is disabled.
%
% 11/26/07    mpr   added timeout if nothing is returned within 10 seconds.
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
% 05/05/17    tbc   Hardcoded PR655 port settings for OSX
% 2018-01-30  tbc   Let FindSerialPort make best guess for port name depending on OS.
%                   (works on OSX & Linux rigs I've used)
%
%

global g_serialPort;
pdur = 0.05; % pause duration between IO 

% This seems to be the default name on OS/X. We don't know about other
% operating systems defaults:
portNumber = FindSerialPort([], 1, 1);

retval = [];
ntries = 0;

% Repeatedly try to establish connection to PR655
while isempty(strfind(retval,'REMOTE')) && ntries<4
    ntries = ntries+1 %#ok<NOPRT>
    
    oldverbo = IOPort('Verbosity', 0);
    pause(pdur);
    
    % usb-serial settings that seem to make PR655 especially happy
    baudRate = 9600;
    parity = 'None';
    dataBits = 8;
    stopBits = 1;
    flowControl = 'None';
    % Compile settings string
    portSettings = sprintf('BaudRate=%i Parity=%s DataBits=%i StopBits=%i FlowControl=%s ',...  %  PollLatency=%2.5f StartBackgroundRead=%i ReadFilterFlags=%i
        baudRate, parity, dataBits, stopBits, flowControl);
    
    if IsOSX
        % Must flush on write, ie., not don't flush on write, at least with PR655
        % on OSX 10.10, as reported in forum message #19808 for more reliable
        % connections:
        try
            g_serialPort = IOPort('OpenSerialPort', portNumber, [portSettings, 'Lenient DontFlushOnWrite=0']);
        catch
            IOPort closeall % ...this is aggressive
            pause(pdur*10)
            g_serialPort = IOPort('OpenSerialPort', portNumber, [portSettings, 'Lenient DontFlushOnWrite=0']);
        end
        
    else
        % On at least Linux (status on Windows is unknown atm.), we must not flush
        % on write - the opposite of OSX behaviour (see forum msg thread #15565):
        try
            g_serialPort = IOPort('OpenSerialPort', portNumber, [portSettings, 'Lenient DontFlushOnWrite=1']);
        catch
            IOPort closeall % ...this is aggressive
            pause(pdur*10)
            g_serialPort = IOPort('OpenSerialPort', portNumber, [portSettings, 'Lenient DontFlushOnWrite=1']);
        end
            
    end
    pause(pdur)
    IOPort('Verbosity', oldverbo);
    pause(pdur)
    
    % Put in Remote Mode --No [CR] after 'PHOTO'
    rm = ['PHOTO', char(13)];
    for i = 1:length(rm)
        IOPort('write', g_serialPort, rm(i));
        pause(.02) % not too short, not too long, juuust right
    end
    
    StartTime = GetSecs;
    retval = [];
    while isempty(retval) && GetSecs-StartTime < 10
        retval = PR655read;
        disp(retval)
        pause(pdur*10)
    end
end

if isempty(retval)
    warning('Could not connect to PR655. Try [un]plugging power & usb cables, and/or using a different usb port...its fickle.')
    keyboard
end

end
