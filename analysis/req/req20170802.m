
sessionList = 'NCPSHB';
sl = get_session_list(sessionList);
nt = numel(sl);

Trials = af(@(t)        MTATrial.validate(t),                     sl    );
stc    = cf(@(t)        t.load('stc'),                            Trials);
xyz    = cf(@(t)        preproc_xyz(t,{'SPLINE_SPINE_HEAD_EQD')), Trials);
fxyz   = cf(@(x)        filter(x.copy(),3,2.5,'low'),             xyz   );
fvxy   = cf(@(x)        vel(x,[],[1,2]),                          fxyz  );
for s=1:nt,             fvxy{s}.data(fvxy{s}.data<1e-3,:) = 1e-3; end
         cf(@(v)        set(v,'data',log10(v.data)),              fvxy  );
ang    = cf(@(t,x)      create(MTADang,t,x),                      Trials,xyz);
fet    = cf(@(t)        fet_bref(t),                              Trials);

% JPDF head height versus head speed
% SET jpdf parameters
ind    = cf(@(s)        [s{'a+w'}],                               stc   );

xAxislabel = 'Head Speed';
edx    = repmat({linspace(0.25,2,100)},[1,nt]);
edy    = repmat({linspace(10,2,180)},  [1,nt]);
out    = cf(@(f,v,i,ex,ey)    hist2([f(i,15),v(i,'head_front')],ex,ey), ...
            fet,fvxy,ind,edx,edy);

figure();
imagesc(edx{1},edy{1},cat)
























