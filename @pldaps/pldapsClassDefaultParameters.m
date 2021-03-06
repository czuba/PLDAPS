function s=pldapsClassDefaultParameters(s)

% % New filename: pldapsClassDefaults.m
% Will try to symlink w/ relative path to updated file
% 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% WARNING:  Outdated file naming
% --------------------------------------------- 
% This file is unused & obsolete.
% 
% pldapsClassDefaults.m  is now used to
% initialize a new experiment session object.
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 

error('Wrong pldaps default setup file...look to @pldaps.pldapsClassDefaults.m for current defaults.')

return

% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 












 if nargin<1
	s=struct;
 end

%s.	.
%s.	behavior.
%s.	behavior.	reward.
 s.	behavior.	reward.	defaultAmount = 0.15;

%s.	datapixx.
 s.	datapixx.	enablePropixxCeilingMount = false;
 s.	datapixx.	enablePropixxRearProjection = false;
 s. datapixx.   rb3d = 0;
 s.	datapixx.	LogOnsetTimestampLevel = 0;
 s.	datapixx.	use = true;
 s.	datapixx.	useAsEyepos = false;
 s.	datapixx.	useForReward = false;

%s.	datapixx.	adc.
 s.	datapixx.	adc.	bufferAddress = [ ];
 s.	datapixx.	adc.	channelGains = 1;
 s.	datapixx.	adc.	channelMapping = 'datapixx.adc.data';
 s.	datapixx.	adc.	channelModes = 0;
 s.	datapixx.	adc.	channelOffsets = 0;
 s.	datapixx.	adc.	channels = [ ];
 s.	datapixx.	adc.	maxSamples = 0;
 s.	datapixx.	adc.	numBufferFrames = 600000;
 s.	datapixx.	adc.	srate = 1000;
 s.	datapixx.	adc.	startDelay = 0;
 s.	datapixx.	adc.	XEyeposChannel = [ ];
 s.	datapixx.	adc.	YEyeposChannel = [ ];

%s.	datapixx.	GetPreciseTime.
 s.	datapixx.	GetPreciseTime.	maxDuration = 0.1; % sec
 s.	datapixx.	GetPreciseTime.	optMinwinThreshold = 0.0003; % sec
 % s.datapixx.	GetPreciseTime.	syncmode = [ ]; % no longer used

%s.	display.
 s.	display.	bgColor = [ 0.5000    0.5000    0.5000 ];
 s.	display.	displayName = 'defaultScreenParameters';
 s.	display.	scrnNum = max(Screen('Screens'));
 s.	display.	useOverlay = 1;
 s.	display.	screenSize = screenSizeSelector(s.display.scrnNum);
 s.	display.	heightcm = 45;
 s.	display.	widthcm = 63;
 s.	display.	viewdist = 57;
 s. display.    ipd = 6.5;
 s. display.    multisample = 0;
 s.	display.	colorclamp = 0;
 s.	display.	forceLinearGamma = false;
 s.	display.	normalizeColor = 1;
 s.	display.	sourceFactorNew = 'GL_SRC_ALPHA';
 s.	display.	destinationFactorNew = 'GL_ONE_MINUS_SRC_ALPHA';
 s.	display.	stereoFlip = [ ];
 s.	display.	stereoMode = 0;
 s. display.    crosstalk = 0;
 s.	display.	switchOverlayCLUTs = false;
 s. display.    useGL = false; % flag for custom 3D rendering features
 s. display.    preOpenScreenFxn = [];
 s. display.    postOpenScreenFxn = [];

%s.	eyelink.
 s.	eyelink.	buffereventlength = 30;
 s.	eyelink.	buffersamplelength = 31;
 s.	eyelink.	calibration_matrix = [ ];
 s.	eyelink.	collectQueue = false;  % No Eyelink queue! This is a timesink, no longer recommended (TBC 2018)
 s.	eyelink.	custom_calibration = false;
 s.	eyelink.	custom_calibrationScale = 0.4; % 0.4 ok for >= 55 cm viewing distance on most (propixx) projection setups
 s.	eyelink.	saveEDF = false;     % better off manually transferring by running this at end of session:      pds.eyelink.fetchEdf(p.trial.session.dir)
 s.	eyelink.	use = true;
 s.	eyelink.	useAsEyepos = true;
 s.	eyelink.	useRawData = false;

%s.	git.
 s.	git.	use = false;
 
%s. keyboard.
 s. keyboard.   devIdx = -1; % PTB default to first keyboard detected

%s.	mouse.
 s.	mouse.	initialCoordinates = [];
 s.	mouse.	use = true;
 s.	mouse.	useAsEyepos = false;
 s.	mouse.	useLocalCoordinates = false;
 s.	mouse.	windowPtr = [ ];

%s.	newEraSyringePump.
 s.	newEraSyringePump.	alarmMode = 1;
 s.	newEraSyringePump.	allowNewDiameter = false;
 s.	newEraSyringePump.	diameter = 38; % 38==60 mL syringe; 29.7==20 mL syringe
 s.	newEraSyringePump.	lowNoiseMode = 0;
 s.	newEraSyringePump.	port = '/dev/cu.usbserial';
 s.	newEraSyringePump.	rate = 60; % ml per minute (...formerly 2900 mL/hr)
 s.	newEraSyringePump.	triggerMode = 'ST';
 s.	newEraSyringePump.	use = false;
 s.	newEraSyringePump.	volumeUnits = 'ML';

%s.	pldaps.
 s.	pldaps.	experimentAfterTrialsFunction = [ ];
 s.	pldaps.	eyeposMovAv = 1;
 s. pldaps. lastBgColor = s.display.bgColor;
 s.	pldaps.	finish = Inf;
 % s.	pldaps.	goodtrial = 0; % This is old; now in p.trial.good. Marking for future deletion. --TBC 2017
 s.	pldaps.	iTrial = 0;
 s.	pldaps.	maxPriority = true;
 s.	pldaps.	maxTrialLength = 30;
 s.	pldaps.	nosave = false;
 s.	pldaps.	pass = false;
 s.	pldaps.	quit = 0;
 s.	pldaps.	trialMasterFunction = 'runModularTrial';
 s.	pldaps.	useFileGUI = false;
 s.	pldaps.	useModularStateFunctions = true;

%s.	pldaps.	dirs.
 s.	pldaps.	dirs.	data = '/Data';
 s. pldaps. dirs.   proot = fileparts(fileparts(mfilename('fullpath'))); % proot == PLDAPS root directory
 s.	pldaps.	dirs.	wavfiles = fullfile( s.pldaps.dirs.proot, 'beepsounds'); %'[proot]/beepsounds';

%s.	pldaps.	draw.
%s.	pldaps.	draw.	cursor.
 s.	pldaps.	draw.	cursor.	use = false;

%s.	pldaps.	draw.	eyepos.
 s.	pldaps.	draw.	eyepos.	use = true;

%s.	pldaps.	draw.	framerate.
 s.	pldaps.	draw.	framerate.	location = [ -30   -10 ];
 s.	pldaps.	draw.	framerate.	nSeconds = 3;
 s.	pldaps.	draw.	framerate.	show = false;
 s.	pldaps.	draw.	framerate.	size = [ 20     5 ];
 s.	pldaps.	draw.	framerate.	use = false;

%s.	pldaps.	draw.	grid.
 s.	pldaps.	draw.	grid.	use = true;

%s.	pldaps.	draw.	photodiode.
 s.	pldaps.	draw.	photodiode.	everyXFrames = 17;
 s.	pldaps.	draw.	photodiode.	location = 1;
 s.	pldaps.	draw.	photodiode.	use = 0;

%s.	pldaps.	pause.
 s.	pldaps.	pause.	preExperiment = true;
 s.	pldaps.	pause.	type = 1;

%s.	pldaps.	save.
 s.	pldaps.	save.	initialParametersMerged = 1;
 s.	pldaps.	save.	mergedData = 0;
 s.	pldaps.	save.	trialTempfiles = 1;
%s.	pldaps.	save.	v73 = 0;

%s.	pldaps.	trialStates.
	tsNeg = -1; tsPos = 1;
 s.	pldaps.	trialStates.        trialSetup = tsNeg;             tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.        trialPrepare = tsNeg;           tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.	frameUpdate = tsPos;            tsPos = tsPos+1;
 s.	pldaps.	trialStates.	framePrepareDrawing = tsPos;    tsPos = tsPos+1;
 s.	pldaps.	trialStates.	frameDraw = tsPos;              tsPos = tsPos+1;
 s.	pldaps.	trialStates.	frameGLDrawLeft = tsPos;        tsPos = tsPos+1;
 s.	pldaps.	trialStates.	frameGLDrawRight = tsPos;       tsPos = tsPos+1;
 s.	pldaps.	trialStates.	frameDrawingFinished = tsPos;   tsPos = tsPos+1;
 s.	pldaps.	trialStates.	frameFlip = tsPos;              tsPos = tsPos+1;
 s. pldaps. trialStates.        trialItiDraw = tsNeg;           tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.        trialCleanUpandSave = tsNeg;    tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.	experimentPreOpenScreen = tsNeg;    tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.	experimentPostOpenScreen = tsNeg;   tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.	experimentAfterTrials = tsNeg;      tsNeg = tsNeg-1;
 s.	pldaps.	trialStates.	experimentCleanUp = tsNeg;          tsNeg = tsNeg-1;
 

%s.	plexon.
%s.	plexon.	spikeserver.
 s.	plexon.	spikeserver.	continuous = false;
 s.	plexon.	spikeserver.	eventsonly = false;
 s.	plexon.	spikeserver.	remoteip = 'xx.xx.xx.xx';
 s.	plexon.	spikeserver.	remoteport = 3333;
 s.	plexon.	spikeserver.	selfip = 'xx.xx.xx.xx';
 s.	plexon.	spikeserver.	selfport = 3332;
 s.	plexon.	spikeserver.	use = 0;
 s. plexon. spikeserver.    spikeCount = 0;

%s.	session.
 s.	session.	experimentSetupFile = [ ];

%s.	sound.
 s.	sound.	deviceid = [ ];
 s.	sound.	use = true;
 s.	sound.	useForReward = true;
 
% Matlab-side [eye]tracker calibration
 s. tracking. use               = true;
 s. tracking. on                = false;
 s. tracking. gridSz            = [20, 16];  % [x,y] size of target grid, in deg
 s. tracking. gridScale         = 1; % scale of calbiration targets relative to display dimensions
 s. tracking. tform             = []; % geometric transform object (index if binocular) See:  pds.tracking.runCalibrationTrial>>updateCalibTransform
 

end



% % % % % % % % % 
% Sub-functions
% % % % % % % % % 

%% screenSizeSelector
%  set default screenRect to 'picture-in-picture' if maxScreens == 0
function screenRect = screenSizeSelector(scrnNum)
if ~scrnNum && scrnNum==max(Screen('Screens'))
    screenRect = floor(Screen('Rect',scrnNum).*0.6)+30;
else
    screenRect = [ ];
end
end
