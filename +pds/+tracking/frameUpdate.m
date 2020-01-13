function p = frameUpdate(p)
% p = pds.tracking.frameUpdate(p)
% get current sample from tracked source & apply calibration
%
% 
% TODO:  Add explainer for how updateFxn should be constructed for different
% tracker sources.
% 
% %!%!%!% WARNING %!%!%!%
% 
%   This method DOES NOT RECORD thefull array of [sample] data
%   - Continuous data recording functions are left up to 
%     tracking source specific code; e.g. pds.eyelink.getQueue
%
% %!%!%!% WARNING %!%!%!%
% 
% See also:  pds.tracking.runCalibrationTrial
% 
% 
% 2020-01-xx  TBC  Wrote it.

% Evolved from pds.eyelink.getQueue.m


if p.trial.tracking.use                     % this is wasteful convention; if you don't use, don't call outside fxn in first place.
    src = p.trial.tracking.source;
    cm = p.trial.(src).calibration_matrix;
    % contents of this should be a geometric transform object,
    % but can also be a simple 3-by-3 transformation matrix (..for now)
    
    % Get current position data
    %     posRaw = pds.mouse.updateFxn(p);
    posRaw = feval(p.static.tracking.updateFxn.(src), p);

    % NOTE: function handles MUST be stored in p.static, NOT in p.trial
    %       ...krufty holdover of the 'params' class; nixing it has long been on the TODO list
    
    % Apply separate calibration matrix to each eye (bino compatible)
    for i = 1:size(posRaw,2)
        if ndims(cm)<3 % isTform
            pos(:,i) = transformPointsInverse(cm(i), posRaw(:,i)')';
            % TODO:  Had to switch to Inverse from transformPointsForward method
            % because forward method doesn't exist for polynomial class tforms
            % ...should check that manual matrix transform is still consistent
            % (or just remove it entirely)
            
        else
            eXY = cm(:,:,i) * [posRaw(1,i), posRaw(2,i), 1]';
            pos(1,i) = eXY(1);
            pos(2,i) = eXY(2);
        end
    end
    
    % % % %Print raw pos to command window while sanity checking
    % % %         fprintf(repmat('\b',1,44));
    % % %         fprintf('%20.18g,  %20.18g\n', posRaw(1:2));
    
    % assign positions to source outputs
    p.trial.tracking.pos = pos;
    p.trial.tracking.posRaw = posRaw;
    %     % Store tracking timecourse w/in source module ??
    %     p.trial.(src).pos(:,p.trial.iFrame) = pos(:);
    %     p.trial.(src).posRaw(:,p.trial.iFrame) = posRaw(:);
    
    % Apply calibrated data as current eye position
    if p.trial.(src).useAsEyepos
        p.trial.eyeX = pos(1,:)';
        p.trial.eyeY = pos(2,:)';
    end
                
end
    
end %main function