
linkSession('Ed10-20140820',...
            '/storage/gravio/data/processed/xyz/Ed10/',...
            '/storage/eduardo/data/processed/nlx/Ed10/');

Session = MTASession('Ed10-20140820',...
                     'cof',...
                      true,...
                     '0x0002',...
                     'vicon',...
                     'nlx',...
                      119.881035...                     
);

Session = MTASession('Ed10-20140820','cof');


QuickTrialSetup(Session,'all',[7,-7]);

Trial = MTATrial('Ed10-20140820');
Trial.stc.updateMode('default');
Trial.stc.load;
Trial.stc.states = {};

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

if isempty(Trial.stc.gsi('h')),
    xyz = Trial.load('xyz');
    xyz.filter('ButFilter',3,2.4);
    fvxy = xyz.vel(3,[1,2]);
    fvxy.data(fvxy.data<1e-3)=1e-3;
    fvxy.data = log10(fvxy.data);
    vper = ThreshCross(fvxy.data,0.5,round(.25*xyz.sampleRate));
    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       vper,...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'velHthresh','h');
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

Trial.stc.save(1);

units = 1:160;
states = {'theta','velthresh','velHthresh'};
overwrite = true;
nsts = numel(states);

for i = 1:nsts,
    pfs{i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[2.2,2.2]);
end

hfig = figure;
unit = units(1);
while unit~=-1,

    for i = 1:nsts,
        subplot2(1,nsts,1,i);cla
        try,
            pfs{i}.plot(unit,[],true,[],false);
            title([pfs{i}.session.trialName ' ' pfs{i}.parameters.states,': ',num2str(unit)]);
        end
    end


    unit = figure_controls(hfig,unit,units);
end



