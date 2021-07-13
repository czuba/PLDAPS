function modNamesOut = updateModNames(p, stimName)
% function pldaps.updateModNames(p, stimName)
% 
% Update list of active and/or particular module names in the pldaps structure.
% - if called w/o output arguements, will assign updated module names to p.trial.pldaps.modNames
%   ELSE, will return modNames struct without making any changes to [p]
% - Optional [stimName] input will constrain matrixModules to those that contain(..., stimName)
% 
% See also: pldaps.getModules, pldaps.run
% 
% 2019-07-26  TBC  Wrote it.
% 2021-07-13  TBC  New matrixModule constraints for compatibility with multiple condMatrix objects
%           

% Establish list of all module names    (see help pldaps.getModules)
modNames.all             = getModules(p, 0);
modNames.matrixModule    = getModules(p, bitset(0,2));
modNames.tracker         = getModules(p, bitset(0,3));

% Constrain matrixModules to ...
if nargin<2
    if isprop(p.condMatrix, 'stimName') && ~isempty(p.condMatrix.stimName)      % 
        % ...the current condMatrix
        stimName = p.condMatrix.stimName;
    elseif ~isfield(p.trial.pldaps.modNames,'currentStim')
        % ...the .currentStim
        stimName = p.trial.pldaps.modNames.currentStim;
    else
        stimName = [];
    end
end


if ~isempty(stimName)
    stimName = cellstr(stimName);
    % match module names to stimName (**if cell of multiple strings, any match will suffice)
    mm = contains(modNames.matrixModule, stimName);
    % remove any unmatched modules
    modNames.matrixModule(~mm) = [];
    % update currentStim with stimName
    modNames.currentStim = stimName;
end

% Combine struct with existing p.trial.pldaps.modNames
% - allowing any additional user-defined module names to carry through
modNames = setstructfields(p.trial.pldaps.modNames, modNames);

% Special case module names get assigned to associated PLDAPS components
if isfield(modNames, 'behavior')
    % extract just this behavior module string (used heavily in behavioral control module [.pmBase])
    p.trial.behavior.modName = modNames.behavior{1};
end

% IF no outputs, assign [modNames] to PLDAPS object
% ELSE just return [modNames] cell
if nargout<1
    % assign to pldaps struct & return
    p.trial.pldaps.modNames = modNames;
    return
else
    % return cell of modNames
    modNamesOut = modNames;
end
    
end %main function

