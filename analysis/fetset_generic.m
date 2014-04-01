function fet = fetset_generic(Trial)

Trial.load('ang');
Trial.load('xyz');

dwin= gausswin(61)./sum(gausswin(61));
fwin = gausswin(11)./sum(gausswin(11));

Trial.xyz.filter(fwin);

fpv = cat(1,[-2,-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','spine_upper','head_front'},[1,2]))));
fpz = cat(1,[-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','head_back'},[3]))));

fet =[];
fet = [fet,fpv,fpz];
fet = [fet,Trial.xyz(:,7,3)-Trial.xyz(:,1,3)];
fet = [fet,Trial.ang(:,3,4,2)];
fet = [fet,Trial.ang(:,5,7,2)];
fet = [fet,fet.^3];
fet = MTADxyz([],[],Filter0(fwin,fet),Trial.xyz.sampleRate);
fet.data(fet.data==0) = nan;


