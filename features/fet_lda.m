function fet = fet_lda(Trial)

xyz = resample(Trial.load('xyz').filter(gtwin(0.75,Trial.xyz.sampleRate)),30);
ang = create(Trial.ang.copy,Trial,xyz);

xyz.addMarker('fhcom',[.7,1,.7],...
              {{'head_back','head_front',[0,0,1]}},...
              xyz.com(xyz.model.rb({'head_left','head_front','head_right'})));

xyz.addMarker('fbcom',[.7,0,.7],...
              {{'spine_lower','pelvis_root',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root'})));


xyz.addMarker('fscom',[.7,0,.7],...
              {{'spine_middle','spine_upper',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_middle','spine_upper'})));

fpv = xyz.vel({'fbcom','fscom','fhcom'},[1,2]);
fpv.data(fpv.data<.01) = .01;

dsv = [0;Filter0(ones(1,7)./7,abs(diff(circ_dist(ang(:,'spine_lower','spine_middle',1),ang(:,'spine_middle','head_front',1)))))];


fet = MTADxyz('data',[fpv.data,...
                      xyz(:,'fhcom',3)-xyz(:,'fbcom',3),...
                      ang(:,'spine_middle','spine_upper',2),...
                      ang(:,'head_back'   ,'head_front' ,2),...
                      ang(:,'spine_lower' ,'spine_upper',3),...
                      ang(:,'spine_lower','spine_middle',3),...
                      ang(:,'pelvis_root','spine_middle',3),...
                      dsv],...
              'sampleRate',xyz.sampleRate);

fet.data(fet.data==0) = nan;