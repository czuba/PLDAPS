function [] = leap(p,state)
    
    switch state
        
        case p.trial.pldaps.trialStates.frameUpdate;
            frameUpdate(p);
                
        case p.trial.pldaps.trialStates.trialCleanUpandSave;
            cleanUpAndSave(p);
        
        case p.trial.pldaps.trialStates.trialSetup;
            trialSetup(p);    
    end

end


function []= trialSetup(p)
    %start listening
    pds.leap.matleap.matleap(2);
    
    %allocate data
    %assuming leap sampling rate o 120Hz for now. Could query later if needed.
    p.trial.leap.cursorSamples = nan(3,round(round(p.trial.pldaps.maxFrames/p.trial.display.frate*120*1.1)));
    p.trial.leap.samplesTimes = nan(3,round(round(p.trial.pldaps.maxFrames/p.trial.display.frate*120*1.1)));
    p.trial.leap.samples = 0;
    
    for iTry=1:10
        if p.trial.leap.samples~=0
            p.trial.leap.succ = true;
            break;
        end
        WaitSecs(0.05);
        getRawCoords(p);
    end
end


function []= frameUpdate(p)
    getRawCoords(p);
end

function [] = cleanUpAndSave(p)
    %stop listening
    pds.leap.matleap.matleap(3);
    %prune unused data
    p.trial.leap.cursorSamples(:,p.trial.leap.samples+1:end) = [];
    p.trial.leap.samplesTimes(:,p.trial.leap.samples+1:end) = [];

end

function getRawCoords(p)
    succ = false;
    f = pds.leap.matleap.matleap(1);
    fpcell={f.pointables};
    if ~isempty(f)
        has_pointables=~cellfun(@isempty, fpcell);
        nsamples=sum(has_pointables);
        if nsamples>0
            succ=true;
            pos=cellfun(@(x) x(1).position, fpcell(has_pointables), 'UniformOutput', false);
            pos=vertcat(pos{:})';
            p.trial.leap.cursorSamples(1:3,p.trial.leap.samples+(1:nsamples))=pos;
            p.trial.leap.samplesTimes(1:3,p.trial.leap.samples+(1:nsamples))=[f(has_pointables).getsecs; f(has_pointables).timestamp; f(has_pointables).id];
            p.trial.leap.samples=p.trial.leap.samples+nsamples;
        end
    end
    p.trial.leap.succ = succ;
end

