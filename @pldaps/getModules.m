function [moduleNames, moduleFunctionHandles, moduleRequestedStates] = getModules(p, moduleType)
% function [moduleNames, moduleFunctionHandles, moduleRequestedStates, moduleLocationInputs] = getModules(p, [moduleType])
% 
% Parse current module activity state, execution order, and trialStates during which they should be active.
% Return a list of module names, handles, states, and 
% This is called at the beginning of each trial.
%
% Module order will be determined by .stateFunction.order, running from -Inf : 0 : Inf.
%       Recommend using "order 0" for the pldapsDefaultTrialFunction module, then ordering other
%       modules before or after it in time.
%       *** NOTE: This is distinct from the ordering of trial states (.stateFunction.requestedStates)
%       becasue negative module order values are not ignored during the trial like state values are.
%       ...this is unfortunate, and potentially confusing.
% 
% [moduleType] is optional parameter for constraining the list of modules returned.
%       def:  moduleType = bitset(0,1); % only active modules
%       Operates as binary flags on each bit in moduleType, where bit:
%       1 == active state of module:        p.trial.(sn).use;           bitget(moduleType,1)
%       2 == marked for use by condMatrix:  p.trial.(sn).matrixModule;  bitget(moduleType,2)
%       3 == tracking module                p.trial.(sn).useAsEyepos;   bitget(moduleType,3)
%            (e.g. .eyelink, or .mouse, ...)
% 
%       EXAMPLES)
%       % Retrieve all active modules (classic/default mode)
%       getModules(p);
% 
%       % Retrieve [tracking module]:
%       getModules(p, bitset(0,3));     % also equivalent to:  getModules(p,3);
% 
%       % Retrieve only [active] [matrixModules]:
%       getModules(p, bitset(bitset(0,1), 2));
% 
%       ....
%       Pro-EXAMPLE) Retrieve modules containing a specific field & value:
%       getModules(p, {'thisField',thisValue; 'thatField',thatValue});
%
% See also:  pldaps.run
% 
% 2016-xx-xx  jk  Written for pldaps 4.2 by Jonas
% 2017-10-xx  tbc Added flag for module active state
% 2018-06-26  tbc Expanded with binary flag moduleType, for more control over module sub-selection & .condMatrix
% 2019-07-23  tbc Updated subselection options
% 2021-07-09  tbc Excised "moduleLocationInputs"...**all** PLDAPS modules must accept standard 3 inputs:  (p, state, sn)
% 

% Parse inputs & Initialize
[moduleNames, moduleFunctionHandles, moduleRequestedStates] = deal([]);
doFieldCheck = 0;

if nargin<2 || isempty(moduleType)
    moduleType = bitset(0,1); % return only 'active' modules:  p.trial.(sn).use == true;
end


%% Find all modules in p.trial struct
moduleNames=fieldnames(p.trial);

% Shortcircuit special cases (allow detection of core PLDAPS fields, like .mouse or .eyelink)
% --- return only 'active' module marked for use as eye position:  p.trial.(sn).use == true & p.trial.(sn).useAsEyepos == true;
if iscell(moduleType)
    % check for presence & value of a specified subfield
    for i = 1:size(moduleType,1)
        thisField = moduleType{i,1};
        % "isfield" call will inherently exclude all non-structs
        moduleNames(cellfun(@(x) ~isfield(p.trial.(x), thisField), moduleNames)) = []; % check for presence
        if size(moduleType,2)>1
            thisVal = moduleType{i,2};
            if ischar(thisVal)
                % use string for comparison to avoid comparison of each char in string
                thisVal = string(thisVal);
            end
            moduleNames = moduleNames(cellfun(@(x) p.trial.(x).(thisField) == thisVal, moduleNames)); % matching value
        end
    end
    return

elseif bitget(moduleType, 3)
    moduleNames = moduleNames(cellfun(@(x) ((isfield(p.trial.(x),'useAsEyepos') && p.trial.(x).useAsEyepos) && (isfield(p.trial.(x),'use') && p.trial.(x).use)), moduleNames));
    return
end

% Characteristics of a PLDAPS module, a field of [p.trial] must
% - be a struct itself
% - contain a [.stateFunction] field 
%   - the ".stateFunction" is [a poorly named/antequated] struct field containing PLDAPS module internals
%   - when the helper function pldapsModule() is used for module creation,
%     most users will never need to see/touch/understand .stateFunction internals
moduleNames(cellfun(@(x) (~isfield(p.trial.(x),'stateFunction')), moduleNames))=[];

% Break free from a zillion calls to the input PLDAPS object [p] by extracting
% the ".stateFunction" struct from all candidate modules
% - slow overhead from excessive [p.] calls due to inefficiencies in PLDAPS "params" hierarchy class
SF = cell2mat(cellfun(@(x) p.trial.(x).stateFunction, moduleNames, 'uni',0));
% Now "arrayfun(...,SF)" replaces old calls to "cellfun(...p.trial.(s).stateFunction, moduleNames...)"

% Subselect modules based on moduleType
% --- return only 'matrixModules'
if bitget(moduleType, 2)
    drop = arrayfun(@(x) (~isfield(x,'matrixModule') || ~x.matrixModule), SF);
    SF(drop) = [];
    moduleNames(drop) = [];
end

% --- return only 'active' modules:  p.trial.(sn).use == true;
if bitget(moduleType, 1)
    drop = cellfun(@(x) (~isfield(p.trial.(x),'use') || ~p.trial.(x).use), moduleNames);
    SF(drop) = [];
    moduleNames(drop) = [];
end


%% Sort module execution by requested .stateFunction.order:
% - module code will be executed in order from -Inf:Inf (...Nan even after Inf, but don't do that)
% - If no .order specified, defaults to Inf.
moduleOrder = inf(size(moduleNames));
% if .order defined, retrieve it
hasOrder = arrayfun(@(x) isfield(x,'order'), SF);

moduleOrder(hasOrder) = arrayfun(@(x) x.order, SF(hasOrder));
[moduleOrder, so] = sort(moduleOrder);
SF = SF(so);
moduleNames = moduleNames(so);


if nargout>1
    % Outside of PLDAPS under-the-hood calls, the rest of this is superfluous
    % (...& full of super slow/klunky Legacy code)
    %
        
    %% Handles to module code
    % - because normal "@" function handles can't exist in PLDAPS params hierarchy, this string to function conversion is a workaround
    moduleFunctionHandles=arrayfun(@(x) str2func(x.name), SF, 'UniformOutput', false);
    
    
    %% Limit module execution to certain trialStates
    % Cross reference all trial states in use, with those specifically requested by the module
    % If none specified, make module active for all states.
    availableStates = fieldnames(p.trial.pldaps.trialStates);
    %a little too long, ok, so if requestedStates is not defined, or .all
    %is true, we will call it for all states. Otherwise it will only call
    %the ones defined and true. --jk 2016(?)
    %

    % Excised mmmany repetitive PLDAPS object usage w/in this crazy call
    requestedStates = cellfun(@(x)...
        (arrayfun(@(y)...
        (~isfield(y,'requestedStates') || strcmpi(y.requestedStates, 'all')...
        || isfield(y.requestedStates,'all')...
        || (isfield(y.requestedStates, x) && y.requestedStates.(x))),... % end @(y) customFxn
        SF)),... % end @(x) customFxn
        availableStates, 'UniformOutput', false);
    
    % -----------------------------------------------------------------------
    % not totally clear what this special case is...backwards compatibility?
    if isfield(p.trial.pldaps,'trialFunction') && ~isempty(p.trial.pldaps.trialFunction)
        keyboard
        % ****************************************************************
        % Encountering this is sign of VERY OUTDATED CODE that should be 
        % retired or updated...bypass this roadblock at your own risk.
        % ****************************************************************
        moduleNames{end+1}='stimulus';
        moduleFunctionHandles{end+1}=str2func(p.trial.pldaps.trialFunction);
        for iState=1:length(requestedStates)
            requestedStates{iState}(end+1)=true;
        end
        moduleOrder(end+1)=NaN; %#ok<NASGU>
    end
    % -----------------------------------------------------------------------
    
    
    % Format requested states to a struct of each available trialState, with logical flags for each module
    % No, logical flags force use of a "find" call every time a module is called on every state (mmany times per frame!!)
    % ...instead, do the find just once here
    requestedStates=cellfun(@(x) find(x), requestedStates, 'UniformOutput', false);
    % NOTE: Something very strange happens here...
    %   Somewhere in the mergeToSingleStruct method of Params class,
    %   1-by-n fields get flipped into n-by-1 (and vice versa!@)
    %   So when this function is called it returns differently shaped inputs
    %   depending on whether its outside or inside the trialMasterFunction  (i.e. experimentPostOpenScreen vs. trialSetup)
    %   ...a proper fix to mergeToSingleStruct should be determined --TBC Oct. 2017
    %
    % The following line normalizes output, but does not address the cause.
    %   (...only trial & error to find out why this line was here in first place...--TBC)
    requestedStates=cellfun(@(x) reshape(x,1,numel(x)), requestedStates, 'UniformOutput', false);
    
    % assign module indices to struct of availableStates
    moduleRequestedStates=cell2struct(requestedStates,availableStates);

end

end %main function
