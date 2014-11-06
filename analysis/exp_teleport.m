%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
%MTAstartup;
%Trial = MTATrial('er06-20130614','fly');
%Session = MTASession('er06-20130614');
%Trial = MTATrial('jg05-20120317');

Trial = MTATrial('Ed10-20140820','all','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telcrtl1','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telshift1','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telcrtl2','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telshift2','rov');
% $$$ Trial = MTATrial('Ed10-20140820','teleport','rov');

tnames = {'telcrtl1','telshift1','telcrtl2','telshift2','teleport'};
nt = numel(tnames);

states = {'theta'};
nsts = size(states,2);

display = true;
overwrite = true;
units = 1:100;

[accg,tbin] = autoccg(Trial,units,'theta');


overwrite = true;
pfs = {};

for t = 1:nt
    Trial = MTATrial('Ed10-20140820',tnames{t},'rov');    
    for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
        %pfs{t,i} = MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
    end
end

if display,

    hfig = figure(666999);
    unit = units(1);
    while unit~=-1,
% $$$         for t = 1:nt,
% $$$         for i = 1:nsts,
% $$$             %subplot2(nt,nsts,t,i);cla
% $$$             subplotfit(t,6);cla
% $$$             %imagesc(pfs{1}.adata.bins{1},pfs{1}.adata.bins{2},reshape(pfs{t}.data.rateMap(:,unit),35,60)),colorbar
% $$$             pfs{t}.plot(unit,[],true);
% $$$             title([pfs{t}.session.trialName ':' pfs{t}.parameters.states,': ',num2str(unit)]);
% $$$         end
% $$$         end
        subplot2(5,2,1,1); cla; pfs{1}.plot(unit,[],true);
        title([pfs{1}.session.trialName ':' pfs{1}.parameters.states,': ',num2str(unit)]);
        subplot2(5,2,2,1); cla; pfs{2}.plot(unit,[],true);
        title([pfs{2}.session.trialName ':' pfs{2}.parameters.states,': ',num2str(unit)]);
        subplot2(5,2,3,1); cla; pfs{3}.plot(unit,[],true);
        title([pfs{3}.session.trialName ':' pfs{3}.parameters.states,': ',num2str(unit)]);
        subplot2(5,2,4,1); cla; pfs{4}.plot(unit,[],true);
        title([pfs{4}.session.trialName ':' pfs{4}.parameters.states,': ',num2str(unit)]);
        subplot2(5,2,5,1); cla; pfs{5}.plot(unit,[],true);
        title([pfs{5}.session.trialName ':' pfs{5}.parameters.states,': ',num2str(unit)]);
        subplot2(5,2,1,2); cla,bar(tbin,accg(:,unit));axis tight;
% $$$         subplotfit(6,6);cla,bar(tbin,accg(:,unit));axis tight;
        reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
                  ['exp_teleport-' Trial.name] ,false,Trial.name);
        unit = figure_controls(hfig,unit,units);
    
    end

end

