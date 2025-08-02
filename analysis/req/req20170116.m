sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
FigDir = create_directory('/storage/gravio/figures/placefields'); 

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};


tind = 6:14;

cf(@(t)  pfs_2d_theta(t,'overwrite',true),  Trials(tind));

units = cf(@(t)  select_placefields(t),  Trials(tind));

cf(@(t,u,s)  pfs_2d_states(t,u,[],s,true,'',true),  Trials(tind),units,repmat({states},[1,numel(tind)]));
