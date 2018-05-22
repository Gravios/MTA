





sessionList = get_session_list('MjgER2016');

Trials = af(@(t) MTATrial.validate(t), sessionList);
         cf(@(t)  t.load('nq'),         Trials    );
units  = cf(@(t) select_placefields(t), Trials    );

units = req20180123_remove_bad_units(units);

sum(cellfun(@numel,newUnits))


% ASSERT that the neuron quality struct and spike cluster map are updated
cellfun(@(t) assert(size(t.nq.eDist,1)==size(t.spk.map,1)), Trials);



tind = 20;
% GENERATE unit autocorrelograms
[accg,tbins] = cf(@(t)  autoccg(t),  Trials(tind));
% LOAD theta state placefields
pft = cf(@(t,u) pfs_2d_theta(t,u,'numIter',1), Trials(tind), units(tind));
% LOAD behavior state placefields
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta','pause&theta','lpause&theta','hpause&theta'};
pfs = cf(@(t,u) pfs_2d_states(t,u,[],[],false,'',false,1), Trials(tind), units(tind));

states = {'loc','lloc','hloc','rear','pause','lpause','hpause'};
pfb = cf(@(t,u,s) pfs_2d_states(t,u,[],s,false,'',false,1), Trials(tind), units(tind),{states});

[pfd,tags] = req20180123_ver5(Trials{tind});

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta','pause&theta','lpause&theta','hpause&theta'};
states = cf(@(t,s) [t.stc{s}]-(t.stc{'r'}+[-1,1]),repmat(Trials(tind),size(states)),states);
states{4} = Trials{tind}.stc{'r'}+[0.5,-0.5];
pfr = cf(@(t,u,s) pfs_2d_states(t,u,[],s,false,'',true,1), Trials(tind), units(tind),{states});
    

nx = numel(pft{1})+numel(pfs{1})+1;


hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position(1:2) = [1, 1];
hfig.Position(3:4) = [nx*5, 9];
hfig.PaperPositionMode = 'auto';

Trial = Trials{tind};
for unit = units{tind},

    clf();
    sp = gobjects([1,0]);
    uicontrol('Parent',  hfig,                                                ...
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
              'Position',[3,1,6,4]);

    % ACCG
    sp(end+1) = subplot(3,nx,1);
    bar(tbins{1},accg{1}(:,unit));
    xlim(tbins{1}([1,end]));
    title(num2str(unit));

    % PFT
    sp(end+1) = subplot(3,nx,2);
    plot(pft{1},unit,1,'text',[],true);
    clear_axes_ticks(sp(end))

    % PFD
    sp(end+1) = subplot(3,nx,2+nx);
    plot(pfd{1},unit,1,'text',[],false);
    clear_axes_ticks(sp(end))
    
    % PFS
    for s = 1:numel(pfs{1}),
        sp(end+1) = subplot(3,nx,2+s);
        plot(pfb{1}{s},unit,1,'text',[],true);
        title(pfr{1}{s}.parameters.states);
        clear_axes_ticks(sp(end))        
        
        sp(end+1) = subplot(3,nx,2+s+nx);
        plot(pfs{1}{s},unit,1,'text',[],true);
        title(pfs{1}{s}.parameters.states);
        clear_axes_ticks(sp(end))        
        
        sp(end+1) = subplot(3,nx,2+s+nx*2);
        plot(pfr{1}{s},unit,1,'text',[],true);
        title(pfr{1}{s}.parameters.states);
        clear_axes_ticks(sp(end))        
    end
    
    af(@(h) set(h,'Units','centimeters'), sp);
    af(@(h) set(h,'Position',[h.Position(1:2),3,3]), sp);

    
    waitforbuttonpress();
end

