
Trials = af(@(s) MTATrial.validate(s), get_session_list_v3('BehaviorPlaceCode'));
tid = 32;
Trial = Trials{tid};
Trial.load('stc','msnn_ppsvd_raux');
pft = pfs_2d_theta(Trial);
pfs = pfs_2d_states(Trial,Trial.spk.map(:,1)');


numStates = numel(pfs);
% SELECT placefields based on Visual Characteristics
su = [];
bu = [];
pu = [];
iu = [];
bh = [];
eg = [];

hfig = figure(7);
for u = Trial.spk.map(:,1)'
    clf(hfig);
    subplot(1,numStates+1,1);
    plot(pft,u,1,'text');
    title(['Unit: ',num2str(u)])
    for sts = 1:numStates,
        subplot(1,numStates+1,sts+1);
        plot(pfs{sts},u,1,'text');
    end
    waitforbuttonpress();
    switch hfig.CurrentCharacter
      case 'e',
        pu(end+1) = u;
        eg(end+1) = u;
      case 'a',
        pu(end+1) = u;
        bh(end+1) = u;
        eg(end+1) = u;
      case 'b',
        pu(end+1) = u;
        bh(end+1) = u;
      case 'p', % good
        pu(end+1) = u;
      case 'i', % good
        iu(end+1) = u;
      case 'x', % bad
        bu(end+1) = u;
      case 's',
        su(end+1) = u;
    end
end

Trial.spk.set_unit_set(Trial,'ego',sort(unique(eg)));
Trial.spk.set_unit_set(Trial,'bhv',sort(unique(bh)));
Trial.spk.set_unit_set(Trial,'placecells',sort(unique(pu)));
Trial.spk.set_unit_set(Trial,'badcells',sort(unique(bu)));
Trial.spk.set_unit_set(Trial,'inteneurons',sort(unique(iu)));
Trial.spk.set_unit_set(Trial,'silentcells',sort(unique(su)));


