
sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
Trials  = af(@(t)  MTATrial.validate(t),             sessionList);
stc     = cf(@(t)  t.load('stc','msnn_ppsvd_raux'),  Trials);


phlr = [];
pdur = [];
for t = 1:numel(Trials),
    pper = stc{t}{'p'};
    plper = cast(stc{t}{'q'},'TimeSeries');
    phper = cast(stc{t}{'j'},'TimeSeries');
    for p = pper.data'
        pls = sum(plper(p'));
        phs = sum(phper(p'));
        phlr(end+1) = pls/(phs+pls);
        pdur(end+1) = diff(p)./pper.sampleRate;
    end
end


xhlr = [];
xdur = [];
for t = 1:numel(Trials),
    xper = stc{t}{'x'};
    xlper = cast(stc{t}{'l'},'TimeSeries');
    xhper = cast(stc{t}{'h'},'TimeSeries');
    for p = xper.data'
        xls = sum(xlper(p'));
        xhs = sum(xhper(p'));
        xhlr(end+1) = xls/(xhs+xls);
        xdur(end+1) = diff(p)./xper.sampleRate;
    end
end


figure();
subplot(221);
hist2([xhlr'+randn([numel(xhlr),1])/100,log10(xdur')],linspace(-0.1,1.1,50),linspace(-2,3,50))
xlabel('low to high ratio');
ylabel('duration (s)');
title('low/high ration vs pause period duration')
caxis([0,50]);
subplot(222);
hist(xhlr'+randn([numel(xhlr),1])/100,100);
xlabel('low to high ratio');
ylabel('count');
title('low/high ratio walk');
subplot(223);
hist2([phlr'+randn([numel(phlr),1])/100,log10(pdur')],linspace(-0.1,1.1,50),linspace(-2,3,50))
xlabel('low to high ratio');
ylabel('duration (s)');
title('low/high ration vs walk period duration')
caxis([0,50]);
subplot(224);
hist(phlr'+randn([numel(phlr),1])/100,100);
title('low/high ratio walk');
xlabel('low to high ratio');
ylabel('count');
