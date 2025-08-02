%cf(@(p) set(p,'tag',p.tag(5:end)), pfd);



cf(@(p) p.update_filename(), pfd(:,[1,2,4,5]));
cf(@(p) p.save()           , pfd(:,[1,2,4,5]));


Trials = af(@(t) MTATrial.validate(t), sessionList);
         cf(@(t)  t.load('nq'),        Trials     );
nunits = cf(@(t) select_placefields(t), Trials);
ounits = units



newUnits = req20180123_remove_bad_units(newUnits);

sum(cellfun(@numel,newUnits))


punits = cf(@(n,u) setdiff(n,u), newUnits,units);

% GENERATE unit autocorrelograms
[accg,tbins] = cf(@(t)  autoccg(t),  Trials);

% LOAD theta state placefields
pft = cf(@(t,u) pfs_2d_theta(t,u,'numIter',1), Trials, punits);
% LOAD behavior state placefields
pfs = cf(@(t,u) pfs_2d_states(t,u,'numIter',1), Trials, punits);

nx = numel(pft{1})+numel(pfs{1})+1+1;



tind = 20;

hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, nx*5, 6];
hfig.PaperPositionMode = 'auto';

Trial = Trials{tind};
for unit = punits{tind},
    clf();
    sp = gobjects([1,0]);
    FigInfo = uicontrol('Parent',  hfig,                                                ...
                        'Style',   'text',                                              ...
                        'String',  {['Unit: ',num2str(unit)],                           ...
                                    Trial.filebase,                                     ...
                                    ['stcMode: ',Trial.stc.mode],                       ...
                                    ['eDist:   ',num2str(Trial.nq.eDist(unit))],        ...
                                    ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],...
                                    ['SNR:     ',num2str(Trial.nq.SNR(unit))],          ...
                                    ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],       ...
                                    ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]  ...
                                   },                                                   ...
                        'Units',   'centimeters',                                       ...
                        'Position',[2,1,6,3.5]);

    % ACCG
    sp(end+1) = subplot(1,nx,2);
    bar(tbins{tind},accg{tind}(:,unit));
    xlim(tbins{tind}([1,end]));
    title(num2str(unit));

    % PFT
    sp(end+1) = subplot(1,nx,3);
    plot(pft{tind},unit,1,text,[],true);
    clear_axes_ticks(sp(end))
    % PFS
    for s = 1:numel(pfs{tind}),
        sp(end+1) = subplot(1,nx,3+s);
        plot(pfs{tind}{s},unit,1,'text',[],true);
        title(pfs{tind}{s}.parameters.states);
        clear_axes_ticks(sp(end))        
    end
    
    af(@(h) set(h,'Units','centimeters'), sp);
    af(@(h) set(h,'Position',[h.Position(1:2),3,3]), sp);

    
    waitforbuttonpress();
end

