% My wife hates me day :(

Trial = MTATrial.validate('jg05-20120310.cof.all');
units = Trial.spk.map(:,1)';
pfs{1} = pfs_2d_states(Trial,units,[],[],[],'debursted',1,1);
pfs{2} = pfs_2d_states(Trial,units,[],[],[],'allspikes',1,1,'all');
pfs{3} = pfs_2d_states(Trial,units,[],[],[],'blge2',1,1,'blge2');
pfs{4} = pfs_2d_states(Trial,units,[],[],[],'blge3',1,1,'blge3');

cf(@(p) p.purge_savefile(), pfs{3})
cf(@(p) p.purge_savefile(), pfs{4})

ny = numel(pfs);
nx = numel(pfs{1});
mrt = max(cell2mat(cf(@(p) max(cell2mat(cf(@(f) f.maxRate(),  p)),[],2),  pfs)),[],2);
mrtw = [0.95,0.95,0.25,0.125];
hfig = figure(39394);
%set(hfig,'Position',[10,30,1560,207])
clf();
for u = units
    for y = 1:ny,
        for x = 1:nx,
            subplot2(ny,nx,y,x);
            plot(pfs{y}{x},u,1,true,mrt(u).*mrtw(y),true);
            title(num2str(u));
        end
    end
    waitforbuttonpress();
end
 
