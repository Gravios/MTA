function gen_unit_profile(Trial,varargin)
[units,states,overwrite,mode,stc_mode,display] = DefaultArgs(varargin,{'pyr',[],false,'pfs','auto_wbhr',0});

if isempty(states),   states = Trial.stc.list_state_attrib('key'); end
numsts = numel(states);

Trial.stc.updateMode(stc_mode);Trial.stc.load;

if isempty(Trial.nq), Trial.load('nq');                     end
if Trial.xyz.isempty, Trial.load('xyz');                    end

if ischar(units),    units = select_units(Trial,18,units); end




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

for s = 1:numsts,
    [accg(:,:,s),tbin] = autoccg(Trial,units,states{s});
end


figH = figure(2324);
for u = units,
    %clims = [0, max(cellfun(@maxRate,pfs,repmat({u},1,numsts)))];
clims = [0,30];
s = 1;
    for s = 1:numsts,

        subplot2(2,numsts,1,s);
        pfs{s}.plot(u,'colorLimits',clims);
        title([pfs{s}.parameters.states ' ' num2str(u)]);

        subplot2(2,numsts,2,s);
        bar(tbin,accg(:,u==units,s)),axis tight,

    end

waitforbuttonpress
savefig(fullfile(Trial.spath,'figures',[Trial.filebase '-unit_prof-' num2str(u)]),figH,'png');
end

