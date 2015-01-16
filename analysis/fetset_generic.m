function fet = fetset_generic(Trial)

Trial = MTATrial('jg05-20120310');
ang = Trial.ang.copy;
xyz = Trial.load('xyz');

dwin = gtwin(.5,xyz.sampleRate);
fwin = gtwin(.1,xyz.sampleRate);

xyz.filter(fwin);
ang.create(Trial,xyz);

vel = xyz.vel({'spine_lower','spine_upper','head_front'},[1,2]);
vel.filter(dwin);
fpv = cat(1,vel.data);

vel = xyz.vel({'spine_lower','head_front'},[3]);
vel.filter(dwin);
fpz = cat(1,vel.data);
rolf = fet_roll(Trial,[],'default');

fet =[];
fet = [fet,fpv,fpz];
fet = [fet,xyz(:,'head_front',3)-xyz(:,'spine_lower',3)];
fet = [fet,ang(:,'pelvis_root','spine_middle',2)];
fet = [fet,ang(:,'spine_middle','spine_upper',2)];
fet = [fet,ang(:,'head_back','head_front',2)];
fet = [fet,circ_dist(ang(:,'spine_middle','spine_upper',1),ang(:,'head_back','head_front',1))];
fet = [fet,rolf.data];

xyz.filter(dwin);
ang.create(Trial,xyz);

fet = [fet,abs([0;diff(ang(:,'spine_lower','spine_middle',2))])];
fet = [fet,abs([0;diff(ang(:,'head_back','head_front',2))])];


nz = nniz(fet);
fet(nz,:) = unity(fet(nz,:));


% $$$ mset = [1,2];
% $$$ xl = -1.5:.03:1.5;
% $$$ yl = -.6:.05:1.4;
% $$$ s = 'w';
% $$$ figure,hist2([rolf(nniz(rolf)),ang(nniz(rolf),mset(1),mset(2),2)],100,100);caxis([0,1000])
% $$$ figure,hist2([xyz(nniz(rolf),1,3),ang(nniz(rolf),mset(1),mset(2),2)],100,100);caxis([0,100])
% $$$ figure,hist2([xyz(nniz(rolf),1,3),ang(nniz(rolf),mset(1),mset(2),2)],xl,yl);caxis([0,10000])
% $$$ figure,hist2([xyz(Trial.stc{s},1,3),ang(Trial.stc{s},mset(1),mset(2),2)],xl,yl);caxis([0,1000])
% $$$ sd = xyz.copy;
% $$$ sd.data = circ_dist(ang(:,'spine_middle','spine_upper',1),ang(:,'head_back','head_front',1));
% $$$ figure,hist2([sd(nniz(rolf)),ang(nniz(rolf),mset(1),mset(2),2)],xl,yl);caxis([0,4000])
% $$$ figure,hist2([sd(Trial.stc{s}),ang(Trial.stc{s},mset(1),mset(2),2)],xl,yl);caxis([0,400])



% $$$ figure,imagesc((1:size(fet,1))/xyz.sampleRate,1:size(fet,2),fet'),caxis([-2,3])
% $$$ Lines(Trial.stc{'r'}(:)/xyz.sampleRate,[],'g',[],5); Lines(Trial.stc{'w'}(:)/xyz.sampleRate,[],'k',[],5);
% $$$ colormap hot
% $$$  
% $$$ fet = MTADxyz('data',fet,'sampleRate',Trial.xyz.sampleRate);
% $$$ 
% $$$ nz = Trial.stc{'h'}.copy;
% $$$ I = [];
% $$$ for i = 1:fet.size(2),
% $$$     for j = 1:fet.size(2),
% $$$         h = hist2([fet(nz,i),fet(nz,j)],-4:.2:4,-4:.2:4);
% $$$         h = h/sum(h(:));
% $$$         lh1 = log(sum(h,1));
% $$$         lh2 = log(sum(h,2));
% $$$         I(i,j) = nansum(nansum(h .* bsxfun(@minus,bsxfun(@minus,log(h),lh1),lh2) ));        
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ [rhm,fs,ts]  = fet_rhm (Trial,[],'wcsd');
% $$$ %figure,imagesc(ts,fs,log10(rhm.data)'),caxis([-6,-4]),axis xy
% $$$ swag = fet_swag(Trial,[],'wcsd');
% $$$ roll = fet_roll(Trial,[],'wcsd');
% $$$ 
% $$$ [U,S,V] = svd(unity(rhm(Trial.stc{'w'},1:5:end)));
% $$$ urhm = unity(rhm(:,1:5:end));
% $$$ for i=1:rhm.size(1),score(i,:) = V(:,1:4)'*urhm(i,:)';end
% $$$ 
% $$$ 
