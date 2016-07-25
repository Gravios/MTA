MTAstartup('vr_exp');

Session = MTASession('Ed10-20140821',...
                     'rov',...
                      true,...
                     '0x0002',...
                     'vicon',...
                     'nlx',...
                     149.9974321...
);

xyz = Session.load('xyz');
xyz.data(:,:,1) = xyz.data(:,:,1)-70;
xyz.data(:,:,2) = xyz.data(:,:,2)-325;
xyz.save;


Trials = SessionList('Ed10VR_opticflow',...
                     '/storage/gravio/data/processed/xyz/Ed10/',...
                     '/storage/eduardo/data/processed/nlx/Ed10/');
%QuickTrialSetup(Trials);

Trial = MTATrial('Ed10-20140821','rov','all');
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
    fvxy = xyz.vel(6,[1,2]);
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


% Calculate and plot
nt = numel(Trials);
Stc = Trial.stc.copy;

states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

display = true;
overwrite = false;
units = 1:160;

[accg,tbin] = autoccg(Trial,units,'theta');

pfs = {};
for t = 1:nt
    Trial = MTATrial(Trials(t).sessionName,...
                     Trials(t).mazeName,...
                     Trials(t).trialName);
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial);
    for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[2.2,2.2]);
    end
end





if display,

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for t = 1:nt,
        for i = 1:nsts,
            subplot2(nt+1,nsts,t,i);cla
            try,pfs{t,i}.plot(unit,[],true,[],false);
            title([pfs{t,i}.session.trialName ' ' pfs{t,i}.parameters.states,': ',num2str(unit)]);
            
            end
        end
        end
        subplot2(nt+1,nsts,t+1,1);cla
        bar(tbin,accg(:,unit));
        %reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
        %          ['exp_optflow-' Trial.name] ,false,Trial.name);

        unit = figure_controls(hfig,unit,units);
    end

end



if display,

    
    spOpts.width  = 4;
    spOpts.height = 2;
    spOpts.ny = numel(tnames)+1;
    spOpts.nx = numel(states);
    spOpts.padding = 2;
    spOpts.units = 'centimeters';

    figOpts.units = 'centimeters';
    figOpts.headerPadding = 4;
    figOpts.footerPadding = 4;
    figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];

    
    sp = [];
    autoincr = true;

    figHnum = 666999;
    set(0,'defaultAxesFontSize',8,...
          'defaultTextFontSize',8)
    hfig = figure(figHnum);clf
    set(hfig,'units',figOpts.units)
    set(hfig,'Position',figOpts.position)
    set(hfig,'PaperPositionMode','auto');

    unit = units(1);
    while unit~=-1,
        clf
        pfsMaxRate = zeros([nt,nsts]);
        for t = 1:nt,
            for i = 1:nsts,                            
                pfsMaxRate(t,i) = pfs{t,i}.maxRate(unit);                
            end
        end
        
        for t = 1:nt,
            for i = 1:nsts,
                sp(t,i) = axes('Units',spOpts.units,...
                               'Position',[(spOpts.width +round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                           (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                                            spOpts.width,...
                                            spOpts.height]...
                );

                pfs{t,i}.plot(unit,[],true,[0,max(pfsMaxRate(:,i))],false);
                title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
            end
        end

        
        t = t+1;
        i = 1;
        sp(t,i) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width +round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                   (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                                    spOpts.width,...
                                    spOpts.height]...
        );
        bar(tbin,accg(:,unit));
        xlim([min(tbin),max(tbin)]);
        title([' AutoCCG: Unit ',num2str(unit)])

        print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'figures/vr_exp/Ed10-20140820-shift_teleport',...
                                     ['pfs_',num2str(unit),'.eps']))
        print(gcf,'-dpng',fullfile(getenv('PROJECT'),'figures/vr_exp/Ed10-20140820-shift_teleport',...
                                   ['pfs_',num2str(unit),'.png']))
% $$$         
% $$$         reportfig('/storage/gravio/figures/',...
% $$$                   hfig,...
% $$$                   ['exp_teleport-' Trial.name],...
% $$$                   'vr_exp',...
% $$$                   false,...
% $$$                   ['Unit: ' num2str(unit)],...  Tag
% $$$                   '',...                Comment
% $$$                   100,...                Resolution
% $$$                   false,...             SaveFig
% $$$                   'png',...             Format
% $$$                   8,...                 Width
% $$$                   12,...                 Height
% $$$                   unit...               Id
% $$$         );
        unit = figure_controls(hfig,unit,units,autoincr);
    
    end

end
