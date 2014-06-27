function gen_unit_profile(Trial,varargin)
[units,states,overwrite,mode,stc_mode,display] = DefaultArgs(varargin,{'pyr',[],true,'pfs','auto_wbhr',0});


stc_mode = 'auto_wbhr';
% $$$ stc_mode = 'manual_tmknsrw';
% $$$ stc_mode = 'qda_segmented';
% $$$ stc_mode = 'qda_ext_seg';
% $$$ stc_mode = 'extended';
% $$$ stc_mode = 'mlda_ext';
% $$$ stc_mode = 'mlda_ext_t1';


Trial.stc.updateMode(stc_mode);Trial.stc.load;
if isempty(states),   states = Trial.stc.list_state_attrib('label'); end
numsts = numel(states);
states = states(cellfun(@isempty,regexpi(states,'theta')));
states = cellfun(@cat,repmat({2},1,numel(states)),states,repmat({'&theta'},1,numel(states)),'uniformoutput',false);
if numsts>numel(states), states{end+1} = 'theta';end


if isempty(Trial.nq), Trial.load('nq');                     end
if Trial.xyz.isempty, Trial.load('xyz');                    end
if ischar(units),    units = select_units(Trial,10,units); end




pfs = {};
for s=1:numsts,
    switch mode,

      case 'knnpfs'
        pfs{s} = MTAAknnpfs(Trial,units,states{s},overwrite,'numIter',1, ...
                     'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);

      case 'pfs'
        pfs{s} = MTAApfs(Trial,units,states{s},overwrite);

    end
end

% ACCG for each State
for s = 1:numsts,
    [accg(:,:,s),tbin] = autoccg(Trial,units,states{s},'normalization','hz');
end


% SCCG for each State periods must be > 0.5;
Bccg = {};
for s = 1:numsts,
    Bccg{s} = gen_bhv_ccg(Trial,states{s},1);
end

 
figH = figure(2324);
for u = units,
    %clims = [0, max(cellfun(@maxRate,pfs,repmat({u},1,numsts)))];
clims = [0,25];
s = 1;
    for s = 1:numsts,

        subplot2(4,numsts,1,s);
        bar(tbin,accg(:,u,s)),axis tight,

        subplot2(4,numsts,2,s);
        pfs{s}.plot(u,'colorLimits',clims);
        title([pfs{s}.parameters.states ' ' num2str(u)]);

% $$$         subplot2(4,numsts,3,s);
% $$$         Bccg{s}.plot(u,1); axis tight,
% $$$ 
% $$$         subplot2(4,numsts,4,s);
% $$$         Bccg{s}.plot(u,2); axis tight,

    end

waitforbuttonpress
reportfig(Trial,'FileName',['unit_profile_' stc_mode],'Comment',num2str(u))
%savefig(fullfile(Trial.spath,'figures','unit_profiles',[Trial.filebase '-unit_prof-' num2str(u) '.png']),figH,'png');
%saveas(figH,fullfile(Trial.spath,'figures','unit_profiles',[Trial.filebase '-unit_prof-' num2str(u) '.pdf']));
end
 

