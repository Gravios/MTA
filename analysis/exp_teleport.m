MTAstartup('vr_exp')

% $$$ Session = MTASession('Ed10-20140820',...
% $$$                      'rov',...
% $$$                       true,...
% $$$                      '0x0002',...
% $$$                      'vicon',...
% $$$                      'nlx',...
% $$$                      149.9974321...
% $$$ );
% $$$ Session = MTASession('Ed10-20140820','rov');
% $$$ xyz = Session.load('xyz');
% $$$ xyz.data(:,:,1) = xyz.data(:,:,1)-70;
% $$$ xyz.data(:,:,2) = xyz.data(:,:,2)-325;
% $$$ xyz.save;

QuickTrialSetup('Ed10VR_teleport');



%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
%MTAstartup;
%Trial = MTATrial('er06-20130614','fly');
%Session = MTASession('er06-20130614');
%Trial = MTATrial('jg05-20120317');

Trial = MTATrial('Ed10-20140820','rov','all');
Trial.stc.updateMode('default');
Trial.stc.load;
% $$$ Trial = MTATrial('Ed10-20140820','telcrtl1','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telshift1','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telcrtl2','rov');
% $$$ Trial = MTATrial('Ed10-20140820','telshift2','rov');
% $$$ Trial = MTATrial('Ed10-20140820','teleport','rov');

tnames = {'tel_C1','tel_S1','tel_C2','tel_S2','tel_A1'};


if isempty(Trial.stc.gsi('v')),
    xyz = Trial.load('xyz');
    xyz.filter('ButFilter',3,2.4);
    fvxy = xyz.vel(1,[1,2]);
    fvxy.data(fvxy.data<1e-3)=1e-3;
    fvxy.data = log10(fvxy.data);
    vper = ThreshCross(fvxy.data,0.5,round(.25*xyz.sampleRate));
    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       vper,...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'velthresh','v');
end

if isempty(Trial.stc.gsi('r')),
    rper = rear(Trial,'com',45);
    Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   rper,...
                   xyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'rear','r');
end


if isempty(Trial.stc.gsi('n')),
    Trial.stc.states{end+1} = Trial.stc{'v'}-(Trial.stc{'r',120}+[-.5,.5]);
    Trial.stc.states{end}.key = 'n';
    Trial.stc.states{end}.label = 'NRvel';    
    Trial.stc.states{end}.updateFilename([Trial.filebase,'.sst.',...
                                          Trial.stc.states{end}.label,'.',...
                                          Trial.stc.states{end}.key,'.mat']);
end


if isempty(Trial.stc.gsi('t')),
    Trial = labelTheta(Trial,[],32);
end


Stc = Trial.stc.copy;
nt = numel(tnames);
states = {'theta','velthresh','rear','NRvel'};
nsts = size(states,2);

display = true;
overwrite = false;
units = 1:160;

[accg,tbin] = autoccg(Trial,units,'theta');

pfs = {};
stc = Trial.stc.copy;


for t = 1:nt
    Trial = MTATrial('Ed10-20140820','rov',tnames{t});    
    Trial.stc = Stc.copy;
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
                %subplot2(nt,nsts+1,t,i);cla
                subplot2(nt+1,nsts,t,i);cla

                %imagesc(pfs{1}.adata.bins{1},pfs{1}.adata.bins{2},reshape(pfs{t}.data.rateMap(:,unit),35,60)),colorbar
                pfs{t,i}.plot(unit,[],true,[],false);
                title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
            end
        end
        subplot2(nt+1,nsts,nt+1,1);
        bar(accg(:,unit));

        reportfig('/storage/gravio/figures/',...
                  hfig,...
                  ['exp_teleport-' Trial.name],...
                  'vr_exp',...
                  false,...
                  ['Unit: ' num2str(unit)],...  Tag
                  '',...                Comment
                  100,...                Resolution
                  false,...             SaveFig
                  'png',...             Format
                  8,...                 Width
                  12,...                 Height
                  unit...               Id
        );
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
    



