
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
lfp = Trial.lfp.copy;
lfp.filename = 'jg05-20120310.lfp';
lfp.create(Trial,61);
lfp.resample(xyz);
phs = lfp.phase;


states = {'theta','rear&theta','loc&theta'};
nsts = numel(states);
ow = 0;

accg = {};
pfs = {};
spk = {};
for i = 1:nsts,
    pfs{i} = MTAApfs(Trial,[],states{i},ow);
    [accg{i},tibn] = autoccg(Trial,[],states{i},320,25,'hz','deburst');
    spk{i} = Trial.spk.copy;
    spk{i}.create(Trial,xyz.sampleRate,states{i},[],'deburst');
end


lfp = Trial.lfp.copy;
lfp.filename = 'jg05-20120310.lfp';
lfp.create(Trial,61);
tpow = fet_spec(Trial,lfp,'mtcsdglong');
tpow.resample(xyz);


aIncr = false;
hfig = figure(838384);
units = 1:size(Trial.spk.map,1);
unit = 1;
while unit~=-1;
    clf
    for i = 1:nsts,
        try
            subplot2(nsts,4,i,1) ;cla; pfs{i}.plot(unit);  title(pfs{i}.parameters.states)
            subplot2(nsts,4,i,2) ;cla; circ_plot(phs(spk{i}(unit)),'hist',[],30,true,true),title(num2str(unit))
            subplot2(nsts,4,i,3) ;cla; bar(tibn,accg{i}(:,unit)),axis tight,title('accg')
            subplot2(nsts,4,i,4) ;cla; %hist(log10(tpow(spk{i}(unit),15)./tpow(spk{i}(unit),1)),0:.03:3),title('dratio'),xlim([0,3])
        catch
            subplot2(nsts,4,i,1);cla
            subplot2(nsts,4,i,2);cla
            subplot2(nsts,4,i,3);cla
            subplot2(nsts,4,i,4);cla
        end
    end
    
    set(hfig,'papertype','A1');
    set(hfig,'paperposition',get(hfig,'paperposition').*[0,0,0,0]+[0,0,10,7]);
    %    saveas(hfig,['C:\Users\justi_000\Documents\figures\pfs_state_spkphs\' Trial.filebase '.pssp-' num2str(unit) '.png'],'png');
    unit = figure_controls(hfig,unit,units,aIncr);
    %saveas(f,['/gpfs01/sirota/bach/homes/gravio/figures/pfs_state_spkphs/' Trial.filebase '.pssp-' num2str(unit) '.png']);
    
end



