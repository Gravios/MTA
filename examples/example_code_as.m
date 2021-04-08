
% $$$ gravio@theta:~$ mkdir /storage/gravio/data/processed/ephys/as01/
% $$$ gravio@theta:~$ cd /storage/gravio/data/processed/ephys/as01/
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01$ ls
% $$$ as01-20191118
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01$ mkdir as01-20191106
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01$ ln -s /storage2/andrey/data/processed/003282/2019-11-06_22-36-07/* /storage/gravio/data/processed/ephys/as01/as01-20191106/
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01$ cd as01-20191106/
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ ls
% $$$ all_channels.events  all.fet.8  all.res.8        Continuous_Data.openephys
% $$$ all.clu.3            all.h5     all.spk.1        external.traj
% $$$ all.clu.3.bak        all.klg.3  all.spk.2        hotspots.traj
% $$$ all.clu.7            all.klg.7  all.spk.3        messages.events
% $$$ all.clu.7.bak        all.lfp    all.spk.4        positions.traj
% $$$ all.dat              all.nrs    all.spk.5        rewards.traj
% $$$ all.fet.1            all.res.1  all.spk.6        settings.xml
% $$$ all.fet.2            all.res.2  all.spk.7        tracking.short.avi
% $$$ all.fet.3            all.res.3  all.spk.8        units
% $$$ all.fet.4            all.res.4  all.xml          VShiftA.mtl
% $$$ all.fet.5            all.res.5  analysis         VShiftA.obj
% $$$ all.fet.6            all.res.6  arena.traj       VShiftA_prod.json
% $$$ all.fet.7            all.res.7  collisions.traj
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ rename 's/all\./as01-20191106\./' all*
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ ls
% $$$ all_channels.events      as01-20191106.klg.3  as01-20191106.spk.6
% $$$ analysis                 as01-20191106.klg.7  as01-20191106.spk.7
% $$$ arena.traj               as01-20191106.lfp    as01-20191106.spk.8
% $$$ as01-20191106.clu.3      as01-20191106.nrs    as01-20191106.xml
% $$$ as01-20191106.clu.3.bak  as01-20191106.res.1  collisions.traj
% $$$ as01-20191106.clu.7      as01-20191106.res.2  Continuous_Data.openephys
% $$$ as01-20191106.clu.7.bak  as01-20191106.res.3  external.traj
% $$$ as01-20191106.dat        as01-20191106.res.4  hotspots.traj
% $$$ as01-20191106.fet.1      as01-20191106.res.5  messages.events
% $$$ as01-20191106.fet.2      as01-20191106.res.6  positions.traj
% $$$ as01-20191106.fet.3      as01-20191106.res.7  rewards.traj
% $$$ as01-20191106.fet.4      as01-20191106.res.8  settings.xml
% $$$ as01-20191106.fet.5      as01-20191106.spk.1  tracking.short.avi
% $$$ as01-20191106.fet.6      as01-20191106.spk.2  units
% $$$ as01-20191106.fet.7      as01-20191106.spk.3  VShiftA.mtl
% $$$ as01-20191106.fet.8      as01-20191106.spk.4  VShiftA.obj
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv positions.traj  as01-20191118.positions.traj 
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv arena.traj  as01-20191118.arena.traj 
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv messages.events as01-20191118.messages.events
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv as01-20191118.messages.events as01-20191106.messages.events
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv as01-20191118.arena.traj as01-20191106.arena.traj 
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv as01-20191118.positions.traj as01-20191106.positions.traj 
% $$$ gravio@theta:/storage/gravio/data/processed/ephys/as01/as01-20191106$ mv as01-20191118.messages.events as01-20191106.messages.events


sessionList = get_session_list('sobolev');
for sid = 1:numel(sessionList);

    link_session(sessionName,Dpaths);

    s = MTASession(sessionList(sid).sessionName,   ...
                   sessionList(sid).mazeName,      ...
                   true,           ...
                   '',      ...
                   {'openephys','sobolev'},      ...
                   sessionList(sid).xyzSampleRate);

    xyz = s.load('xyz');
    mazeCenter = [sessionList(sid).xOffSet,sessionList(sid).yOffSet];
    rot = sessionList(sid).rotation;
    zoffset = sessionList(sid).zOffSet;
    xyz.data(nniz(xyz),:,:) = bsxfun(@minus,xyz.data(nniz(xyz),:,:),permute([mazeCenter,zoffset],[1,3,2]));
    xyz.data(nniz(xyz),:,[1,2])  = multiprod(xyz.data(nniz(xyz),:,[1,2]),...
                                             [cos(-rot), sin(-rot); -sin(-rot), cos(-rot)],...
                                             [3],[1,2]);
    xyz.save();

    %figure,plot(xyz(:,1,1),xyz(:,1,2),'.')

    NeuronQuality(s);

    % CREATE Trial from Session
    Trial = MTATrial(s.name,'vrr','all',true,xyz.sync.data);

end



for sid = 1:numel(sessionList);

    Trial = MTATrial.validate(sessionList(sid));
    % LABEL theta periods
    labelTheta(Trial,Trial.stc,true,false);


    arena = Trial.load('arena');
    xyz = Trial.load('xyz');



    % CREATE behavioral and maze states
    vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),'hcom',[1,2]);

    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(vxy.data>5,0.5,10),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'vel','v');
    Trial.stc{'v'}.save(1);



    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(xyz(:,3) > 125,0.5,10),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'rear','r');
    Trial.stc{'r'}.save(1);

    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(xyz(:,3) < 57,0.5,10),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'low','l');
    Trial.stc{'l'}.save(1);

    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(xyz(:,3) > 57 & xyz(:,3) < 125,0.5,1),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'high','h');
    Trial.stc{'h'}.save(1);


    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(vxy.data>5 & arena(:,1,1) > -170,0.5,10),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'aori','b');
    Trial.stc{'b'}.save(1);


    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       ThreshCross(vxy.data>5 & arena(:,1,1) < -170,0.5,10),...
                       xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'asft','z');
    Trial.stc{'z'}.save(1);


    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       [Trial.stc{'asft&vel&theta'}.data],...
                       Trial.lfp.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'shift','s');
    Trial.stc{'s'}.save(1);



    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       [Trial.stc{'aori&vel&theta'}.data],...
                       Trial.lfp.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'cntrl','o');
    Trial.stc{'o'}.save(1);


end





pfsSft = {};
pfsOri = {};
% SET placefield args
pfsArgs.binDims          = [25,25];
pfsArgs.SmoothingWeights = [3,3];
pfsArgs.spkMode          = '';
overwrite = true;
for sid = 1:numel(sessionList);
% COMPUTE original maze position placefields
    pfsOri{sid} = MTAApfs(Trial,                                         ... % Trial
                     [],                                               ... % units
                     'cntrl',                                          ... % state 
                     overwrite,                                        ... % overwrite
                     '',                                               ... % tag
                     pfsArgs.binDims,                                  ... % binDims
                     pfsArgs.SmoothingWeights,                         ... % SmoothingWeights
                     'spkMode', pfsArgs.spkMode                        ... % spike mode
                     );
% COMPUTE shifted maze position placefields
    pfsSft{sid} = MTAApfs(Trial,                                         ... % Trial
                     [],                                               ... % units
                     'shift',                                          ... % state 
                     overwrite,                                        ... % overwrite
                     '',                                               ... % tag
                     pfsArgs.binDims,                                  ... % binDims
                     pfsArgs.SmoothingWeights,                         ... % SmoothingWeights
                     'spkMode', pfsArgs.spkMode                        ... % spike mode
                     );
end





% COMPUTE placefield stats
ds = batch_compute_pfstats('sobolev',                           ... % session list name
                           'default',                           ... % state collection name
                           {'cntrl','shift'},                   ... % states
                           'teststats',                         ... % tag 
                           overwrite,                           ... % overwrite
                          'pfsArgs',pfsArgs                     ... % Placefield arguments
);





% PLOT placefields
% keybindings
% q: quit
% i: input unit index
% n: next
% p: previous
hfig = figure();
u = 1;
while u~=-1
    %subplot(6,6,u);
    subplot(121);
    hold('on');    
    pfsOri{sid}.plot(u,1,'text',[],false);
    rectangle('Position',[0-400,180-800,800,1600],'EdgeColor','m','LineWidth',2);
    plot(sq(ds{sid}.pfkstats(1,u).patchCOM(1,1,:,1)),...
         sq(ds{sid}.pfkstats(1,u).patchCOM(1,1,:,2)),'^m');
    subplot(122);    
    hold('on');
    pfsSft{sid}.plot(u,1,'text',[],false);
    rectangle('Position',[0-400,-130-800,800,1600],'EdgeColor','m','LineWidth',2);    
    plot(sq(ds{sid}.pfkstats(2,u).patchCOM(1,1,:,1)),...
         sq(ds{sid}.pfkstats(2,u).patchCOM(1,1,:,2)),'^m');
    title(num2str(u));
    u = figure_controls(hfig,u,Trial.spk.map(:,1));
end
