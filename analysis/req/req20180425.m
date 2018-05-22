

sessionList = get_session_list('MjgER2016');
Trials = af(@(s) MTATrial.validate(s), sessionList(end-2:end))


spkmap_jg05-201203_15-16;
spkmap_jg05-201203_16-17;

pft = cf(@(t)  pfs_2d_theta(t,[],[],[],1), Trials);

figure
k = 0;
n = 18;
for cluGrp = 6:8;
    for u = unit{1}{cluGrp}',
        ou = find(u(2)==unit{2}{cluGrp}(:,1));
        if isempty(ou),continue;end

        mr = maxRate(pft{1},u(1)-1);
        if mr <20,mr = 8;else,mr=30;end
        subplot(n,3,k*3+1);plot(pft{1},u(1)-1,1,false,mr,true);clear_axes_ticks(gca);
        if k == 0, title(Trials{1}.filebase);end
        subplot(n,3,k*3+2);plot(pft{2},u(2),1,false,mr,true);clear_axes_ticks(gca);
        if k == 0, title(Trials{2}.filebase);end        
        if ~isempty(ou),
            subplot(n,3,k*3+3);plot(pft{3},unit{2}{cluGrp}(ou,2),1,false,mr,true);clear_axes_ticks(gca);
            if k == 0, title(Trials{3}.filebase);end
        end
        k = k+1;        
    end
end

