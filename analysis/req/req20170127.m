% req20170127.m 
% segment locomotive trajectories into in/out bound from homebase locations


Trial = MTATrial('jg05-20120309');

state = Trial.stc{'loc'};

[location,occpancy] = identify_homebase_locations(Trial);

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');

xyz.addMarker('homebase',...     Name
              [1,0,0],...  Color
              {{'homebase','hcom',[0,0,255]}},... 
              permute(repmat([location,0],[size(xyz,1),1]),[1,3,2])...
);



ang = create(MTADang,Trial,xyz);

angularOffsetHeadxHombase = MTADxyz('data',circ_dist(ang(:,'acom','bcom',1),ang(:,'acom','homebase',1)),'sampleRate',xyz.sampleRate);


% $$$ figure,plot(angularOffsetHeadxHombase)
% $$$ figure,plot(ang(:,'acom','homebase',3))

nind = nniz(ang(:,'acom','homebase',3));
fdist = MTADxyz('data',nan([size(ang,1),1]),'sampleRate',xyz.sampleRate);
fdist.data(nind) = [diff(ButFilter(ang(nind,'acom','homebase',3),3,1/ang.sampleRate.*0.5,'low'));0];

% $$$ hold on,
% $$$ plot(fdist.data);
% $$$ 
% $$$ figure,
% $$$ plot(diff(fdist());

homebase = state.copy;
homebase.clear;
homebase.label = 'locToHome';
homebase.key   = 'i';
homebase.updateFilename(Trial);

explore = state.copy;
explore.clear;
explore.label = 'locFromHome';
explore.key   = 'o'
explore.updateFilename(Trial);

for period = state.data',
    if sum(fdist(period'))>0,
        explore.data(end+1,:) = period';
    else
        homebase.data(end+1,:) = period';
    end
end


Trial.stc.states{end+1} = homebase;
Trial.stc.states{end+1} = explore;

figure,
eds = linspace(-pi,pi,100);
hm = bar(eds,histc(ang(homebase,'head_back','head_front',2),eds),'histc');
hm.FaceColor = 'c';
hm.FaceAlpha = 0.4;
hm.EdgeColor = 'c';
hm.EdgeAlpha = 0.4;
hold on,
hb = bar(eds,histc(ang(explore,'head_back','head_front',2),eds),'histc');
hb.FaceColor = 'r';
hb.FaceAlpha = 0.4;
hb.EdgeColor = 'r';
hb.EdgeAlpha = 0.4;

pfs_2d_states(Trial,[],{'loc','lloc','locToHome','hloc','locFromHome','rear','pause'});

states = {'loc','rear','pause'};
Spk = {};
Bccg = {};
for s = 1:numel(states),
    Spk{s} = Trial.spk.copy;
    Spk{s}.create(Trial,xyz.sampleRate,[states{s},'&theta'],[],'deburst');
    Bccg{s} = gen_bhv_ccg(Trial,states{s},1);
end


eds = -pi:0.1:pi;
figure,
for unit = 1:110    
    for s = 1:numel(states)
        subplot2(4,3,1,s);cla;
        bar(eds,histc(angularOffsetHeadxHombase(Spk{s}(unit)),eds),'histc');
        subplot2(4,3,2,s);cla;
        bar(eds,histc(angularOffsetHeadxHombase(Trial.stc{states{s}}),eds),'histc');
        subplot2(4,3,3,s);cla;
        Bccg{s}.plot(unit,1);
        subplot2(4,3,4,s);cla;
        Bccg{s}.plot(unit,1);        
    end
    pause(0.4)
end



% First Column: 
%               meanSpikeWaveform, 
%               FullAutoCCG, 
%               text(unit neuronquality stuff),
%               plot(SpkWidthR,AmpSym) 
% First Row  StateAutoCCG
%


set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)


gUnits = 'centimeters';
hfig = figure(20170129)
hfig.PaperPositionMode = 'auto';
hfig.Units = gUnits;
hfig.Position = [0,0,10,10];

tb = annotation('textbox','Units',gUnits,'Position',[2,6,3,3]);
tb.String = {'This is a line','This is the second line'};

hax = axes('Units','centimeters','Position',[6,1,3,3]);

