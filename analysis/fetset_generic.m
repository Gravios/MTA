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


fet =[];
fet = [fet,fpv,fpz];
fet = [fet,xyz(:,'head_front',3)-xyz(:,'spine_lower',3)];
fet = [fet,ang(:,'spine_middle','spine_upper',2)];
fet = [fet,ang(:,'head_back','head_front',2)];
fet = [fet,ang(:,'spine_middle','spine_upper',1)];


% $$$ xyz.filter(dwin);
% $$$ ang.create(Trial,xyz);

fet = [fet,abs([0;Filter0(gausswin(61)./sum(gausswin(61)),diff(ang(:,'spine_lower','spine_middle',2)))])];
fet = [fet,abs([0;Filter0(gausswin(61)./sum(gausswin(61)),diff(ang(:,'head_back','head_front',2)))])];

%fet = [fet,fet.^3];
%fet = MTADxyz([],[],Filter0(fwin,fet),Trial.xyz.sampleRate);
% $$$ fet.data(fet.data==0) = nan;
nz = nniz(fet);
fet(nz,:) = unity(fet(nz,:));



figure,imagesc((1:size(fet,1))/xyz.sampleRate,1:size(fet,2),fet'),caxis([-2,3])
Lines(Trial.stc{'r'}(:)/xyz.sampleRate,[],'g',[],5); Lines(Trial.stc{'w'}(:)/xyz.sampleRate,[],'k',[],5);
colormap hot

fet = MTADxyz('data',fet,'sampleRate',Trial.xyz.sampleRate);

nz = Trial.stc{'w'}.copy;
I = [];
for i = 1:fet.size(2),
    for j = 1:fet.size(2),
        h = hist2([fet(nz,i),fet(nz,j)],-3:.2:3,-3:.2:3);
        h = h/sum(h(:));
        lh1 = log(sum(h,1));
        lh2 = l og(sum(h,2));
        I(i,j) = nansum(nansum(h .* bsxfun(@minus,bsxfun(@minus,log(h),lh1),lh2) ));        
    end
end


rhm  = fet_rhm (Trial,[],'wcsd');
swag = fet_swag(Trial,[],'wcsd');

[U,S,V] = svd(unity(rhm(Trial.stc{'w'},1:5:end)));
urhm = unity(rhm(:,1:5:end));
for i=1:rhm.size(1),score(i,:) = V(:,1:4)'*urhm(i,:)';end


