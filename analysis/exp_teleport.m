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

xyz = Trial.load('xyz');
xyz.filter(gtwin(.5,xyz.sampleRate));
vel = [0;sqrt(sum(diff(xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'}))).^2,3).*xyz.sampleRate./10)];
vel(~nniz(xyz))=0;
vper = ThreshCross(vel,2,round(.5*xyz.sampleRate));


si = Trial.stc.gsi('v');
if isempty(si),
    Trial.stc.states{end+1} = MTADepoch(Trial.spath,Trial.filebase,vper,xyz.sampleRate,Trial.sync.copy,Trial.sync.data(1),'vel','v');
else
    Trial.stc.states{si}.data = vper;
end


rper = rear(Trial,'com',45);

si = Trial.stc.gsi('r');
if isempty(si),
    Trial.stc.states{end+1} = MTADepoch(Trial.spath,Trial.filebase,rper,120,Trial.sync.copy,Trial.sync.data(1),'rear','r');
else
    Trial.stc.states{si}.data = rper;
end

si = Trial.stc.gsi('n');
if isempty(si),
    Trial.stc.states{end+1} = Trial.stc{'v',120}-(Trial.stc{'r',120}+[-.5,.5]);
    Trial.stc.states{end}.key = 'n';
    Trial.stc.states{end}.label = 'NRvel';    
    Trial.stc.states{end}.updateFilename([Trial.filebase,'.sst.',Trial.stc.states{end}.label,'.',Trial.stc.states{end}.key,'.mat']);
else
    Trial.stc.states{si} = Trial.stc{'v',120}-(Trial.stc{'r',120}+[-.5,.5]);
    Trial.stc.states{si}.key = 'n';
    Trial.stc.states{si}.label = 'NRvel';    
    Trial.stc.states{si}.updateFilename([Trial.filebase,'.sst.',Trial.stc.states{end}.label,'.',Trial.stc.states{end}.key,'.mat']);
    Trial.stc.states{si}.data([1,end],:) = [];
end


nt = numel(tnames);
states = {'theta','vel','rear','NRvel'};
nsts = size(states,2);

display = true;
overwrite = false;
units = 1:275;

[accg,tbin] = autoccg(Trial,units,'theta');

pfs = {};
stc = Trial.stc.copy;


for t = 1:nt
    Trial = MTATrial('Ed10-20140820',tnames{t},'rov');    
    Trial.stc = stc.copy;
    Trial.stc.load(Trial);
    for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
        %pfs{t,i} = MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
    end
end

autoincr = true;
if display,

    hfig = figure(666999);
    unit = units(1);
    while unit~=-1,
        for t = 1:nt,
            for i = 1:nsts,
                subplot2(nt,nsts+1,t,i);cla

                %imagesc(pfs{1}.adata.bins{1},pfs{1}.adata.bins{2},reshape(pfs{t}.data.rateMap(:,unit),35,60)),colorbar
                pfs{t,i}.plot(unit,[],true);
                title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
            end
        end
        subplot2(nt,nsts+1,1,nsts+1); cla,bar(tbin,accg(:,unit));axis tight;
% $$$         subplot2(5,2,1,1); cla; pfs{1}.plot(unit,[],true);
% $$$         title([pfs{1}.session.trialName ':' pfs{1}.parameters.states,': ',num2str(unit)]);
% $$$         subplot2(5,2,2,1); cla; pfs{2}.plot(unit,[],true);
% $$$         title([pfs{2}.session.trialName ':' pfs{2}.parameters.states,': ',num2str(unit)]);
% $$$         subplot2(5,2,3,1); cla; pfs{3}.plot(unit,[],true);
% $$$         title([pfs{3}.session.trialName ':' pfs{3}.parameters.states,': ',num2str(unit)]);
% $$$         subplot2(5,2,4,1); cla; pfs{4}.plot(unit,[],true);
% $$$         title([pfs{4}.session.trialName ':' pfs{4}.parameters.states,': ',num2str(unit)]);
% $$$         subplot2(5,2,5,1); cla; pfs{5}.plot(unit,[],true);
% $$$         title([pfs{5}.session.trialName ':' pfs{5}.parameters.states,': ',num2str(unit)]);
% $$$         subplot2(5,2,1,2); cla,bar(tbin,accg(:,unit));axis tight;
% $$$         subplotfit(6,6);cla,bar(tbin,accg(:,unit));axis tight;

        reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
                  ['exp_teleport-' Trial.name] ,false,Trial.name);
        unit = figure_controls(hfig,unit,units,autoincr);
    
    end

end



unts = [4,25,27,234,235];
hfig = figure(666989);
unit = units(1);

i = 4;

for t=1:5,
for unit = unts,
    sp = subplot2(nt,numel(unts),t,find(ismember(unts,unit)));cla
    pfs{t,i}.plot(unit,[],true);
    title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
end
end


reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
          ['exp_teleport_selected_units-' Trial.name] ,false,Trial.name);
    



