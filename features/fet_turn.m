function turn_feature = fet_turn(Trial)
%Trial = MTATrial('jg05-20120317');


fwin = gtwin(1,Trial.xyz.sampleRate);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(fwin);

rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
Trial.addMarker(xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(fwin,hcom),[1,3,2]));

rbb = xyz.model.rb({'spine_lower','pelvis_root'});
bcom = xyz.com(rbb);
Trial.addMarker(xyz,'fbcom',[.7,0,.7],{{'spine_lower','pelvis_root',[0,0,1]}},permute(Filter0(fwin,bcom),[1,3,2]));

rbu = xyz.model.rb({'spine_middle','spine_upper'});
scom = xyz.com(rbu);
Trial.addMarker(xyz,'fscom',[.7,0,.7],{{'spine_middle','spine_upper',[0,0,1]}},permute(Filter0(fwin,scom),[1,3,2]));


ang = Trial.ang.copy;
ang.create(Trial,xyz);


hb_dir_swing = [];
hs_dir_swing = [];
for i= 2:2:(round(Trial.xyz.sampleRate/2)+(mod(round(Trial.xyz.sampleRate/2),2)~=0)), 
    hb_dir_swing(:,i/2) = circshift(abs(circ_dist(ang(:,'fhcom','fbcom',1),circshift(ang(:,'fhcom','fbcom',1),i))),-i/2);
    hs_dir_swing(:,i/2) = circshift(abs(circ_dist(ang(:,'fhcom','fscom',1),circshift(ang(:,'fhcom','fscom',1),i))),-i/2);
end
%figure, imagesc(unity(hb_dir_swing)')

nhbds = Filter0(gtwin(0.1,ang.sampleRate),mean(hb_dir_swing,2).*mean(diff(hb_dir_swing,1,2),2));
dnhbds = abs([0;diff(nhbds)]);
dnhbds(dnhbds<0.0006) = 0.0006;
dnhbds = dnhbds*1e3;
% $$$ dkn = MTADxyz('data',nhbds/.dnhbds,'sampleRate',Trial.xyz.sampleRate);
turn_feature = nhbds./dnhbds;


nhsds = Filter0(gtwin(0.1,ang.sampleRate),mean(hs_dir_swing,2).*mean(diff(hs_dir_swing,1,2),2));
dnhsds = abs([0;diff(nhsds)]);
dnhsds(dnhsds<0.0006) = 0.0006;
dnhsds = dnhsds*1e3;
% $$$ dkn = MTADxyz('data',nhbds/.dnhbds,'sampleRate',Trial.xyz.sampleRate);
turn_s_feature = nhsds./dnhsds;

tbfet = turn_feature;
tsfet = turn_s_feature;

nshift = round(0.35*Trial.ang.sampleRate);
nshift = nshift+(mod(nshift,2)~=0); % make nshift even

etbfet = GetSegs(tbfet,1:size(tbfet),nshift,0); 
etsfet = GetSegs(tsfet,1:size(tsfet),nshift,0); 

ctccg = zeros([nshift,size(tsfet)]);
for s = 1:nshift,
     ctccg(s,:) = mean(etbfet+circshift(etsfet,[1,s-nshift/2]).*abs((etbfet-circshift(etsfet,[1,s-nshift/2]))./(etbfet+circshift(etsfet,[1,s-nshift/2]))));
end

turn_feature = log10(abs(Filter0(fwin,min(ctccg)')));
% $$$ figure,plot(dkn(Trial.stc{'w'}),vel(Trial.stc{'w'}),'.')
% $$$ hold on,plot(dkn(Trial.stc{'n'}),vel(Trial.stc{'n'}),'.g')
