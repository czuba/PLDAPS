classdef pldaps < handle
%pldaps    main class for PLDAPS (Plexon Datapixx PsychToolbox) version 4.2
% The pldaps contructor accepts the following inputs, all are optional, but may be required later for experiments to run:
%     1. a subject identifier (string)
%     2. a function name or handle that sets up all experiement parameters
%     3. a struct with changes to the defaultParameters, this is usefull for debugging, but could also be used to replace the function.
% As long as the inputs are uniquely identifiable, the order is not important, otherwise 
% for the remaining unclassified inputs the above order is assumed.
% Read README.md for a more detailed explanation of the default usage

%% --- Properties ---
 properties
    defaultParameters params

    conditions cell %cell array with a struct like defaultParameters that only hold condition specific changes or additions

    condMatrix % machinery for selection & randomization of parameters set by the contents of p.conditions
    
    trial       %will get all variables from defaultParameters + correct conditions cell merged. This will get saved automatically. 
                %You can add calculated paraneters to this struct, e.g. the
                %actual eyeposition used for caculating the frame, etc.
    data cell
    
    static      % Storeage for handles and other data that needs to remain static across trials.
                % WARNING: contents of .static are outside of the params class (p.trial), therefore
                % changes may not be fully stored/tracked across trials. 
                % Replaces disused .functionHandles from PLDAPS 4.1
 end

 
%% --- Methods ---
 methods
     function p = pldaps(varargin)
         %% setup default parameters
        %classdefaults: create default structure from function
        defaults{1} = pldaps.pldapsClassDefaults;
        defaultsNames{1} = 'pldapsClassDefaults';
        
        %rigdefaults: load from prefs?
        defaults{2} = getpref('pldaps');
        defaultsNames{2} = 'pldapsRigPrefs';
        
        p.defaultParameters=params(defaults,defaultsNames);
        
        %DEPRECATED:   Use createRigPrefs.m to update rig prefs framework
        if isField(p.defaultParameters,'pldaps.rigParameters')
            error(['Storing rigPrefs within the .pldaps.rigParameters field is depreciated.\n',...
                    '\tRun createRigPrefs.m to create updated preferences storage inline with PLDAPS ver. 4.2 (and beyond)'], [])
        end
        
        
        %% Process inputs
        %if an input is a struct, this is added to the defaultParameters. 
        %if an input is a cell. this is set as the conditions
        
        %It's contents will override previous parameters
        %the first nonStruct is expected to be the subject's name
        %the second nonStruct is expected to be the experiment functionname
        structIndex=cellfun(@isstruct,varargin);
        if any(structIndex)
            if sum(structIndex)>1
                error('pldaps:pldaps', 'Only one struct allowed as input.');
            end
            constructorStruct=varargin{structIndex};
        else
            constructorStruct=struct;
        end

        cellIndex=cellfun(@iscell,varargin);
        if any(cellIndex)
            if sum(cellIndex)>1
                error('pldaps:pldaps', 'Only one cell allowed as input.');
            end
            p.conditions=varargin{cellIndex};
        end
        
        if nargin>4
            error('pldaps:pldaps', 'Only four inputs allowed for now: subject, experimentSetupFile (String or function handle), a struct of parameters and a cell with a struct of parameters for each trial.');
        end
        subjectSet=false;
        for i=1:nargin
            if isa(varargin{i}, 'function_handle') %fucntion handle will be the experimentSetupFunction
                 constructorStruct.session.experimentSetupFile=func2str(varargin{i});
            elseif isa(varargin{i}, 'string')
                    constructorStruct.session.subject=varargin{i};
                    subjectSet=true;                
            elseif isa(varargin{i}, 'char')
                if ~subjectSet  %set experiment file
                    constructorStruct.session.subject=varargin{i};
                    subjectSet=true;
                else
                    constructorStruct.session.experimentSetupFile=varargin{i};
                end
            end
        end
        constructorStruct.session.caller = dbstack(1, '-completenames');
        if size(constructorStruct.session.caller, 1)==0
            % outputs of dbstack get weird if empty (and uninterpretable 'params class' error occurs...srsly not going there again)
            constructorStruct.session.caller = '';
        end
        p.defaultParameters.addLevels({constructorStruct, struct},{'ConstructorInputDefaults', 'SessionParameters'});
        
        
        % Establish p.trial as a handle to p.defaultParameters        
        p.trial = p.defaultParameters; 
        % Hackish duplication, but standard procedure evolved to basically only use/interact with p.trial.
        % Explicity use of p.defaultParameters only really done in legacy or under-the-hood code now.
        % Also allows the same code that works inside a running session (or inside a module) to be run
        % in the command window [...for the most part].
        
        % Take module inventory
        if p.trial.pldaps.useModularStateFunctions
            % Establish list of all module names
            p.trial.pldaps.modNames.all = getModules(p, 0);
            p.trial.pldaps.modNames.matrixModule = getModules(p, bitset(0,2));
        end


     end   
     
     
     %% PDS = save(p, savedFileName)
     % Moving output struct operations into PLDAPS class methods

     function PDS = save(p, savedFileName)
         % parse inputs
         if nargin<2 || isempty(savedFileName)
             savedFileName = fullfile(p.trial.session.dir, 'pds', p.trial.session.file);
         end
         
         % create output struct
         PDS = struct;
         
         fn = fieldnames(p);
         % check for presence of a params class object
         for i = 1:length(fn)
             hasParamsClass = isa(p.(fn{i}),'params');
             if hasParamsClass
                 break
             end
         end
         
         for i = 1:length(fn)
             switch class(p.(fn{i}))
                 case 'params'
                     % get the raw contents of Params hierarchy (...not for mere mortals)
                     [rawParamsStruct, rawParamsNames] = p.(fn{i}).getAllStructs();
                     
                     % Partition baseline parameters present at the onset of all trials (*)
                     PDS.pdsCore.initialParameters       = rawParamsStruct(p.static.pldaps.baseParamsLevels);
                     PDS.pdsCore.initialParameterNames   = rawParamsNames(p.static.pldaps.baseParamsLevels);
                     PDS.pdsCore.initialParameterIndices = p.static.pldaps.baseParamsLevels;
                     
                     % Include a less user-hostile copy of  [p.trial]  in output struct
                     % ! ! ! NOTE: baseParamsLevels can be changed during experiment (i.e. during a pause),
                     % ! ! ! so this merged struct could be misleading.
                     % ! ! ! Truly activeLevels are documented on every trial in:  data{}.pldaps.activeLevels
                     oldLevels = p.(fn{i}).setLevels(p.static.pldaps.baseParamsLevels);
                     PDS.baseParams = mergeToSingleStruct(p.(fn{i}));
                     % put things back the way you found them
                     p.(fn{i}).setLevels(oldLevels);
                     
                     levelsCondition = 1:length(rawParamsStruct);
                     levelsCondition(ismember(levelsCondition, p.static.pldaps.baseParamsLevels)) = [];
                     PDS.conditions = rawParamsStruct(levelsCondition);
                     PDS.conditionNames = rawParamsNames(levelsCondition);
                     
                 case 'condMatrix'
                     if ~isempty(p.(fn{i}))
                         PDS.(fn{i}) = p.(fn{i});
                         % clear out GUI handles
                         PDS.(fn{i}).H = [];
                     end
                     
                 case 'struct'
                     switch fn{i}
                         case 'trial'
                             if hasParamsClass
                                 % p.trial output is already converted to PDS.baseParams by params class save
                             else
                                 % treat as any other struct
                                 PDS.(fn{i}) = p.(fn{i});
                             end
                         otherwise
                             PDS.(fn{i}) = p.(fn{i});
                     end
                     
                 otherwise
                     switch fn{i}
                         case 'conditions'
                             if hasParamsClass
                                 % do nothing. ...another krufty params class workaround
                             else
                                 % treat as any other field
                                 PDS.(fn{i}) = p.(fn{i});
                             end
                         otherwise
                             PDS.(fn{i}) = p.(fn{i});
                     end
             end %class switch
         end %p fieldnames
         
         
         fprintLineBreak
         if nargout==0
             % save to file
             save(savedFileName, '-mat', '-struct', 'PDS');
             fprintf('\tPLDAPS data file saved as:\n\t\t%s\n', savedFileName);
         else
             % return struct & save instructions
             saveCmdString = sprintf('save(''%s'', ''-mat'', ''-struct'', ''PDS'')', savedFileName);
             fprintf(2, '\tPLDAPS data output returned as struct (but not saved!)\n\tSave manually:\n\t\t%s\n', saveCmdString);
         end
         fprintLineBreak
             
     end %save
     
     
 end
 

 %% --- Static Methods ---
methods(Static)
    % Status of these conversion functions is unknown, and all display a warning directing users tooward
    % pds.<method> versions instead. For now, I will disable them, & listen for complaints. --TBC Summer 2018
    % % %     [xy,z] = deg2px(p,xy,z,zIsR)
    % % %
    % % %     [xy,z] = deg2world(p,xy,z,zIsR)
    % % %
    % % %     [xy,z] = px2deg(p,xy,z)
    % % %
    % % %     [xy,z] = world2deg(p,xy,z)
    
    [xy,z] = deg2world(p, varargin);%(p,xy,z)
    
    [held, dist] = checkFixation(varargin);
    
    % Shorthand reward adjustments
    moreReward(varargin) % default up by 10%
    lessReward(varargin) % default down by 10%
    setReward(varargin) % default set to 0.15
    % Call from command window with:
    % >> p.moreReward
    
    s = pldapsClassDefaults(s); %Parameters(s)
    
    [stateValue, stateName] = getReorderedFrameStates(trialStates,moduleRequestedStates)
    
end %static methods


end % classdef
