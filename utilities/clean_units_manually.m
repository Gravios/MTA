

Bfs      = cf(@(T,U)  compute_bhv_ratemaps(T,U),          Trials, Units);
mask_bhv =load('/storage/gravio/data/project/general/analysis/pfsHB_mask.mat');

units = Trial.spk.get_unit_set(Trial,'bhv');

bfs = compute_bhv_ratemaps(Trial,Trial.spk.map(:,1)');
tid = 32;
Trial = Trials{tid};
Trials{tid}.load('nq');
nq = Trials{tid}.nq;


[mccg,tbin] = autoccg(Trial, [], 'theta-sit-groom', 10, 64, 'hz');

filter_edist = nq.eDist(units) < 15;

bd = [];
hfig = figure(14);
for unit = units(:)'
    gcf();
    subplot2(2,7,1,1);
        bar(tbin,mccg(:,unit));
    subplot2(2,7,1,2);
        plot( nq.AvSpk(unit,:));
    subplot2(2,7,1,3);
    bfs.plot(unit,1,'colorbar',[],true,[],false,[],@jet,mask_bhv.mask);
    for sts = 1:numStates
        subplot2(2,numStates,2,sts);
        pfs{sts}.plot(unit, 1, 'text','colorMap',@jet);
    end
    waitforbuttonpress();
    switch hfig.CurrentCharacter
      case 'b',
        bd(end+1) = unit;
    end
end

bad_units = Trial.spk.get_unit_set(Trial,'badcells');
bad_units = unique(sort([bad_units, bd, units(filter_edist)]));
Trial.spk.set_unit_set(Trial,'badcells',bad_units);