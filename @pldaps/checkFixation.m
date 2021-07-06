function [held, dists] = checkFixation(p, sn, xyPix)
% function held = checkFixation(p, sn, xyPix)
% 
% Default fixation check modular PLDAPS
% Uses pdist to compare pixel coordinates from:
%    [p.trial.eyeX, p.trial.eyeY]
% with
%     p.trial.(sn).fixPosPx &
%     p.trial.(sn).fixLimPx
% 
% Limits evaluated based on p.trial.(sn).mode:
%   (mode 2 recommended)
%   0 = skip/passthrough, always report held=true
%   1 = square window, fixLimPx=[x,y] or [x]
%       - xy limits are half-width of full fixation window
%       - singular value will give square window
%       - only operates on one eye position (not binocular)
%   2 = circle window, fixLimPx=[x,y] or [x]
%       - xy limits are radius along horizontal & vertical axis
%       - singular value will give circle
%       - operates on all eye positions available (bino or mono)
%         , held=true only if BOTH eyes w/in limits
%   
% Assigns logical fixation held status to:  p.trial.(sn).isheld
% and returns logical [held] variable, if requested
% 
% If no [sn] input, tries to use p.trial.pldaps.modNames.currentFix{1},
% otherwise defaults to 'fix'.
% 
%
% 2018-08-16  TBC  Wrote it.
% 2020-10-13  TBC  Cleaned & commented

% assume true to allow code to run w/o any eye tracking setup
held = true;

if p.trial.mouse.useAsEyepos || p.trial.eyelink.useAsEyepos
    if nargin<2 || isempty(sn)
        try
            sn = p.trial.pldaps.modNames.currentFix{1};
        catch
            sn = 'fix';
        end
    end
    
    if nargin<3
        % use .eyeX & .eyeY as comparison pixel position(s) 
        eyeXY = [ (p.trial.eyeX-p.trial.display.ctr(1)), -(p.trial.eyeY-p.trial.display.ctr(2)) ];
        limMode = p.trial.(sn).mode;
    else
        eyeXY = xyPix;
        limMode = 2; % always use circular limit window if [xyPix] input supplied
    end
    
    switch limMode
        case 2
            % default to circle window limits as [radius, x/y ratio]
            limRat = p.trial.(sn).fixLimPx(:);
            if numel(limRat)==1
                limRat = [1 1];
            else
                limRat = limRat./limRat(1);
            end
            
            fixXY = p.trial.(sn).fixPosPx(1:2);
            % calc eyeXY relative to center of screen
            dists = pdist2( fixXY, eyeXY, 'seuclidean', limRat);
            held = all(dists ...
                   <=  p.trial.(sn).fixLimPx(1));
            
            
            % fprintf('%s\t%s\t%s\n', pad(mat2str(fixXY,4),18), pad(mat2str(eyeXY,4),18), pad(mat2str(dists,4),18));

               
        case 0
            p.trial.(sn).isheld = held;
            return
        otherwise
            % square window  (***this is jank, and not recommended)
            % calc eyeXY relative to center of screen
            eyeXY = [ (p.trial.eyeX-p.trial.display.ctr(1)), -(p.trial.eyeY-p.trial.display.ctr(2)) ];
            held = squarewindow(0, eyeXY - p.trial.(sn).fixPosPx, p.trial.(sn).fixLimPx(1), p.trial.(sn).fixLimPx(end));
            dists = [];
    end
    if nargin<3
        % Record in pldaps structure before returning
        % - only if no alternate input pixel positions supplied
        p.trial.(sn).isheld = held;
    end
end

end %main function
