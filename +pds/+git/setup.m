function p = setup(p)
%function p = pds.git.setup(p)
%
% Create all git tracking status & info necessary to fully reconstitute source code used
% during the experiment
% - enabled by default & called by [pldaps.run.m] during main initialization
% 
% Each fieldname under [p.trial.git] should correspond to a separate git repo to be tracked
% - ...naturally, the [.use] flag will automatically be excluded (enabled by PLDAPS default; .use == 1)
% - [.pldaps]  repo will be tracked by default
% - [.huklabbasics]  is also tracked by default, if not part of your codebase, consider as example
% 
% - additional repos can be added by including repo name & main function names to p.trial.git:
%       p.trial.git.<myRepo>.mainFxn = 'aUniqueFilename.m';
%       e.g.  p.trial.git.pldaps.mainFxn = 'pldaps.m';
% - for additional info, see comments in code
% 
% ---
% [REVIVAL EXAMPLE]
% To revive a git repo to the source/state of a PLDAPS experiment session:
% 
%     % load a saved PLDAPS data file [.PDS]
%     thisFile = fullfile( myDataPath, 'myDataFile.PDS');
%     pds = load(thisFile,'-mat');
% 
%     % create a unique destination for the new repo
%     myNewRepo = fullfile(pwd, sprintf('pldaps_%s', pds.baseParams.session.file(1:end-4)));
% 
%     % clone a copy from the repo origin
%     cloneRepoStr = sprintf('git clone %s %s', pds.baseParams.git.pldaps.remote.origin, myNewRepo);
%     [err, resp] = system( cloneRepoStr )
% 
%     % checkout the appropriate commit revision
%     checkoutStr = sprintf('git -C %s checkout %s', myNewRepo, pds.baseParams.git.pldaps.revision);
%     [err, resp] = system( checkoutStr )
% 
%     % create a patch file from the git diff at the time of PLDAPS execution
%     if ~isempty(pds.baseParams.git.pldaps.diff)
%         fid = fopen('mydiff.patch','w+');
%         fwrite(fid, pds.baseParams.git.pldaps.diff),
%         fwrite(fid, sprintf('\n')), % ensure trailing newline character
%         fclose(fid);
%         % apply patch to your new repo
%         [err, resp] = system(sprintf('git -C %s apply --whitespace=nowarn %s', myNewRepo, '../mydiff.patch'))
%     end
% 
%     fprintf('Done! New repo located in:\n\t%s\n', myNewRepo);
% ---
% 
% 05/2014 jk wrote it
% 2020/04/06  TBC  Modernized and expanded with complete tracking & revival functionality
% 

if ~isField(p.trial,'git.use') || ~p.trial.git.use
    return
end

% Always track Pldaps git
p.trial.git.pldaps.mainFxn = 'pldaps.m';
p.trial.git.pldaps.basePath = p.trial.pldaps.dirs.proot;
% ****************************************************************
% NOTE:  Best practices to NOT hard code [git.<repo>.basePath]
% - PLDAPS repo is unique because we inherently know more about it
% For all other repos:
% - define  [p.trial.git.<yourRepo>.mainFxn]  as the name of a core/unique function in your repo
%   - [.mainFxn] should preferably be in the lowest/base directory of your repo
%   - if it is inside an "@" matlab class (like pldaps.m), this function will try to automatically unwrap it
%   - if [.mainFxn] does not exist or isempty, will attempt to use the repo name itself:   what(<myRepo>)
% - use of  which(mainFxn)  will ensure results of this git tracking operation reflect
%   **the actual code** that exists in your path at time of execution
% ****************************************************************

% find other git repos that should be tracked
fn = fieldnames(p.trial.git);
fn = ['pldaps', fn(~contains(fn, {'pldaps','use'}))];

% Attempt git tracking on each repo
for i = 1:length(fn)
    try
        % determine [basePath] of repo
        if ~isfield(p.trial.git.(fn{i}), 'basePath') || isempty(p.trial.git.(fn{i}).basePath)
            if ~isfield(p.trial.git.(fn{i}), 'mainFxn') || isempty(p.trial.git.(fn{i}).mainFxn)
                % if no main function provided, assume repo name matches a directory in matlab path'
                p.trial.git.(fn{i}).mainFxn = [];
                w = what(fn{i});
                if ~isempty(w.path)
                    p.trial.git.(fn{i}).basePath = w(1).path;
                else
                    p.trial.git.(fn{i}).basePath = []; % still not found
                end
            else
                % find current version in path
                p.trial.git.(fn{i}).basePath = fileparts( which(p.trial.git.(fn{i}).mainFxn) );
            end
            
            % unwrap any enclosing matlab "@" class in the .basePath
            while contains(p.trial.git.(fn{i}).basePath, '@')
                p.trial.git.(fn{i}).basePath = fileparts(p.trial.git.(fn{i}).basePath);
            end
        end
        thisBase = p.trial.git.(fn{i}).basePath; % shorthand
        
        % Retrieve info about repo source, status, version/commit, & diff
        % - [pds.git.git.m] is just a lightweight wrapper for command line git
        p.trial.git.(fn{i}).status          = pds.git.git(['-C ' thisBase ' status']);
        p.trial.git.(fn{i}).remote.origin   = pds.git.git(['-C ' thisBase ' remote get-url origin']);
        p.trial.git.(fn{i}).branch          = pds.git.git(['-C ' thisBase ' symbolic-ref --short HEAD']);
        p.trial.git.(fn{i}).revision        = pds.git.git(['-C ' thisBase ' rev-parse HEAD']);
        p.trial.git.(fn{i}).diff            = pds.git.git(['-C ' thisBase ' diff']);
    catch
        % separate error field allows everything up to error to carry through
        errString = sprintf('%s repo not found, or inaccessible.',fn{i});
        warning(errString);
        p.trial.git.(fn{i}).error      = errString;
    end
end

end %main function
