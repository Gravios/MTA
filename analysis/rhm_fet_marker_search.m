%Trial = MTATrial('Ed05-20140528');
Trial = MTATrial('Ed10-20140814');
%Trial = MTATrial('jg05-20120317');
%stc_mode = 'qda_filtf1p5';
stc_mode = 'auto_wbhr';
Trial.stc.updateMode(stc_mode);Trial.stc.load;

xyz = Trial.xyz.copy;
xyz.load(Trial);


%% Calculate the Center of Mass of the head for each frame
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
% not filtered
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},hcom);
% filtered
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,255]}},ButFilter(hcom,3,[2]./(Trial.ang.sampleRate/2),'low'));



%% Calculate Rotation Matrix: normal of the head plane as axis of rotation


xyz_hb_b = sq(xyz(:,'head_back',:)-xyz(:,'hcom',:));
xyz_hb_r = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
xyz_hb_l = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
xyz_hf_r = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
xyz_hf_l = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );


head_norm = cross(xyz_hb_b,xyz_hb_r);
head_norm = multiprod(head_norm,1./sqrt(sum(head_norm.^2,2)),2);
j =1:3;
head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
j = [ 0,-1, 1;...
      1, 0,-1;...
     -1, 1, 0];
k = [1,3,2;...
     3,1,1;...
     2,1,1];
head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
rot_ang = deg2rad(45);
head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])+sin(rot_ang)*head_cpm+(1-cos(rot_ang))*head_kron;

j =1:3;

% Rotated marker;
nmark = permute(sum(head_rotMat.*permute(reshape(xyz_hb_b(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2]);

xyz.addMarker('head_br45',[.7,0,.7],{{'head_back','head_right',[0,0,255]}},permute(nmark,[1,3,2])+hcom)



plot3(xyz_hb_b(1:1000,1),xyz_hb_b(1:1000,2),xyz_hb_b(1:1000,3)),
old on,
plot3(xyz_hb_r(1:1000,1),xyz_hb_r(1:1000,2),xyz_hb_r(1:1000,3),'r'),
plot3(nmark(1:1000,1),nmark(1:1000,2),nmark(1:1000,3),'g')



txyz = xyz.data;
txyz = txyz-repmat(txyz(:,10,:),[1,size(txyz,2),1]);
cxbr1 = (txyz(:,8,:)+txyz(:,5,:))./repmat(sqrt(sum((txyz(:,8,:)+txyz(:,5,:)).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3])+xyz(:,10,:);
cxbr1(~nniz(cxbr1(:)))=0;
xyz.addMarker(['hm_' num2str(1)],[.7,0,.7],{{'head_back','head_front',[0,0,1]}},cxbr1);

cxbr2 = (txyz(:,8,:)+txyz(:,7,:))./repmat(sqrt(sum((txyz(:,8,:)+txyz(:,7,:)).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3])+xyz(:,10,:);
cxbr2(~nniz(cxbr2(:)))=0;
xyz.addMarker(['hm_' num2str(2)],[.7,0,.7],{{'head_back','head_front',[0,0,1]}},cxbr2);

txyz = xyz.data;
txyz = txyz-repmat(txyz(:,10,:),[1,size(txyz,2),1]);
cxbr3 = (txyz(:,13,:)+txyz(:,7,:))./repmat(sqrt(sum((txyz(:,13,:)+txyz(:,7,:)).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3])+xyz(:,10,:);
cxbr3(~nniz(cxbr3(:)))=0;
xyz.addMarker(['hm_' num2str(3)],[.7,0,.7],{{'head_back','head_front',[0,0,1]}},cxbr3);

txyz = xyz.data;
txyz = txyz-repmat(txyz(:,10,:),[1,size(txyz,2),1]);
ixbr1 = cross(txyz(:,5,:),txyz(:,8,:),3)./repmat(sqrt(sum(cross(txyz(:,5,:),txyz(:,8,:),3).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3]);
cxbr4 = (ixbr1+txyz(:,14,:))./repmat(sqrt(sum((ixbr1+txyz(:,14,:)).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3])+xyz(:,10,:);
cxbr4(~nniz(cxbr4(:)))=0;
xyz.addMarker(['hm_' num2str(3)],[.7,0,.7],{{'head_back','head_front',[0,0,1]}},cxbr4);






%pis = [-pi:8:pi,-pi:8:pi,0];
% $$$ pis = sin(circ_dist(sq(ang(:,'hcom','head_back',:)),sq(ang(:,'hcom','head_right',:))));
% $$$ pis = pis./repmat([4,4,0],[size(pis,1),1]);
% $$$ pis(~nniz(pis(:))) = 0;

% $$$ for i = 1:8,
% $$$ apis = asin(pis*i);
% $$$ apis(:,3)=0);
% $$$ vMar = sq(ang(:,'hcom','head_back',:))-apis;
% $$$ [X,Y,Z] = sph2cart(vMar(:,1),vMar(:,2),vMar(:,3));
% $$$ vMar = cat(3,X,Y,Z)+xyz(:,'hcom',:);
% $$$ h = line([xyz.data(21,[10],1),vMar(21,1,1)]',[xyz.data(21,[10],2),vMar(21,1,2)]',[xyz.data(21,[10],3),vMar(21,1,3)]'),set(h,'color',[1,0,0])
% $$$ %xyz.addMarker(['hm_' num2str(i)],[.7,0,.7],{{'head_back','head_front',[0,0,1]}},vMar);
% $$$ end

% $$$ txyz = xyz.data;
% $$$ txyz = txyz-repmat(xyz(:,10,:),[1,11,1]);
% $$$ brx = (txyz(:,8,:)+txyz(:,5,:))./repmat(sqrt(sum((txyz(:,8,:)+txyz(:,5,:)).^2,3)),[1,1,3]).*repmat(sqrt(sum(txyz(:,5,:).^2,3)),[1,1,3])+xyz(:,10,:);;

tt = 1040;
figure,
h = line(xyz.data(tt,[5,10],1),xyz.data(tt,[5,10],2),xyz.data(tt,[5,10],3)),set(h,'color',[0,1,0])
line(xyz.data(tt,[6,10],1),xyz.data(tt,[6,10],2),xyz.data(tt,[6,10],3))
line(xyz.data(tt,[7,10],1),xyz.data(tt,[7,10],2),xyz.data(tt,[7,10],3))
line(xyz.data(tt,[8,10],1),xyz.data(tt,[8,10],2),xyz.data(tt,[8,10],3))
%h = line([xyz.data(tt,[10],1),vMar(tt,1,1)]',[xyz.data(tt,[10],2),vMar(tt,1,2)]',[xyz.data(tt,[10],3),vMar(tt,1,3)]'),set(h,'color',[1,0,0])
%h = line([xyz.data(tt,[10],1),brx(tt,1,1)]',[xyz.data(tt,[10],2),brx(tt,1,2)]',[xyz.data(tt,[10],3),brx(tt,1,3)]'),set(h,'color',[1,0,0])

for i = 1:5;
h = line(xyz.data(tt,[10,11+i],1),xyz.data(tt,[10,11+i],2),xyz.data(tt,[10,11+i],3)),set(h,'color',[1,0,0])
end


nang = Trial.ang.copy;
nang.create(Trial,xyz);


%rhm = fet_rhm(Trial);
%rhm = [bang,lang,rang,fang];

% $$$ for i = 1:4,
% $$$ rhm(:,i) = ButFilter(nang(:,['hm_',num2str(i)],'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
% $$$ end
rhm = [];
rhm(:,1) = ButFilter(nang(:,'head_back','fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
rhm(:,2) = ButFilter(nang(:,12,'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
rhm(:,3) = ButFilter(nang(:,8,'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
rhm(:,4) = ButFilter(nang(:,13,'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
rhm(:,5) = ButFilter(nang(:,14,'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
rhm(:,6) = ButFilter(nang(:,15,'fhcom',3),3,[1,30]./(nang.sampleRate/2),'bandpass');
%ncp = fet_ncp(Trial,'chans',[1,2]);
wang = [rhm,ncp];
wang = WhitenSignal(wang,[],1);

[ys,fs,ts] = mtcsdglong(wang,2^9,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);


szy = size(ys);
padding = zeros([round(2^6/xyz.sampleRate/diff(ts(1:2))),szy(2:end)]);
ys = MTADlfp('data',cat(1,padding,ys,padding),'sampleRate',1/diff(ts(1:2)));

% Speed of the head
vxyz = Trial.xyz.copy;
vxyz.load(Trial);
vxyz.filter(gtwin(1,Trial.xyz.sampleRate));
vh = vxyz.vel('spine_lower',[1,2]);

%% Bug: zeros should remain after resample, anti-alias filter seems
%% to be affecting the edges
vh.resample(ys);
vh.data(~nniz(vh)) = nan;


%edges_labels = mat2cell(edges,2,ones(1,size(edges,2)))';
chan_labels = {'head back','head left','head right','head front','Nasal Epithelium Signal','Nasal Cavity Pressure'};

%for s = 1:numel(Trial.stc.states);
s = 'w';

sind = Trial.stc{s};
%sind = Trial.stc.states{s};

ind = nniz(vh(sind));
vhs = abs(log10(vh(ind)));
vhlim =prctile(vhs,[5,98]);
edges = linspace(vhlim(1),vhlim(2),9);
edges = [edges(1:end-1);edges(2:end)];
edges_labels = mat2cell(10.^mean(edges),1,ones(1,size(edges,2)))';
yss = ys(sind,:,:,:);
yss = yss(ind,:,:,:);
clear ind;

vsc = [];
for i = edges,
    ind = i(1)>=vhs&vhs<i(2);% & 5>=mean(log10(yss(:,fs>13&fs>5,3,3)),2);
    for j = 1:size(yss,3),
        for k = 1:size(yss,3),
            if sum(ind)~=0,
                if k~=j,
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))./mean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
                else
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(yss(ind,:,k,j));
                end

            else
                vsc(find(i(1)==edges(1,:)),:,:) = zeros([1,size(yss,2),size(yss,3)]);
            end

        end
    end
end


% $$$ for j = 1:size(yss,3),
% $$$     for k = 1:size(yss,3),
rhm_mar = [5,12,8],%13,14,15];
k = size(yss,3);
figure,
for j = 1:size(yss,3)-2,

tt = 1040;

        subplot2(6,size(yss,3)-2,[1,2],j),
h = line(xyz.data(tt,[5,10],1),xyz.data(tt,[5,10],2),xyz.data(tt,[5,10],3)),set(h,'color',[0,1,0])
line(xyz.data(tt,[6,10],1),xyz.data(tt,[6,10],2),xyz.data(tt,[6,10],3))
line(xyz.data(tt,[7,10],1),xyz.data(tt,[7,10],2),xyz.data(tt,[7,10],3))
line(xyz.data(tt,[8,10],1),xyz.data(tt,[8,10],2),xyz.data(tt,[8,10],3))
h = line(xyz.data(tt,[rhm_mar(j),10],1),xyz.data(tt,[rhm_mar(j),10],2),xyz.data(tt,[rhm_mar(j),10],3)),set(h,'color',[1,0,0])
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])



        subplot2(6,size(yss,3)-2,3,j),
        hist(vhs,linspace(vhlim(1),vhlim(2),9)),axis tight,
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])

        subplot2(6,size(yss,3)-2,[4:6],j),
        if j~=k
            imagesc(1:size(edges,2),fs,vsc(:,:,j,k)'),axis xy,
            if j~=1,set(gca,'yticklabel',[]),end
            %title([chan_labels{j} ' & ' chan_labels{k} ' coherence'])
        else
            imagesc(1:size(edges,2),fs,log10(vsc(:,:,j,k))'),axis xy,
            %title([chan_labels{j} ' PSD X Body Speed'])
        end
        set(gca,'XTickLabel',cellfun(@num2str,cellfun(@round,cellfun(@transpose,edges_labels(2:2:8),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)');
caxis([.3,.8])
        %if j~=k,caxis([.3,.9]),end
        %if j==size(yss,3),xlabel('Binned Body Speed cm/s'),end
        %if k==1,ylabel('Frequency Hz'),end
        %colorbar
% $$$     end
end
%suptitle(['Coherence During ' sind.label]);

colorbar

%ForAllSubplots('colorbar')


reportfig(Trial, ...
          'FileName','mean_Coherence_RHM_NCP','Comment',['State: ' ...
                    Trial.stc.states{s}.label]  )

end




% phase difference vs mean ncp[6,12] power

mar = 5;%'hm_1';
myncp = median(log10(ys(:,12>fs&fs>6,end,end)),2);
myrhm = median(log10(ys(:,12>fs&fs>6,mar,mar)),2);
mynrp = circ_median(angle(ys(:,12>fs&fs>6,mar,end)),[],2);
sper = Trial.stc{'r',ys.sampleRate}.cast('TimeSeries');
pind = nniz(myrhm)&nniz(myncp)&sper.data(1:end-1);

%ysl = 3.3;
%sind = sper.data(1:end-1)&mys>ysl;
%figure,histcirc(circ_mean(angle(ys(sind,12>fs&fs>6,mar,end)),[],2),30)


plims = prctile(myncp(pind),[2,98]);
pedgs = linspace(2,plims(2),10);
[~,pbins] = histc(myncp(pind),pedgs);

phdgs = linspace(-pi,pi,30);
[~,hbins] = histc(mynrp(pind),phdgs);

A = accumarray([pbins+1,hbins+1],1,[10,30]);
%figure,
imagesc(phdgs(2:end),pedgs(2:end),A(2:end,2:end)),axis xy
ylabel('ncp 6-12Hz pow');
xlabel('phase difference (ncp-rhm)');
title(['ncp-rhm, power X phase diff, JPDF: ' sper.label])


figure,hist2([clip(myncp(pind),0,8),clip(myrhm(pind),-8,-3)],40,40)
caxis([0,150])
xlabel('ncp 6-12Hz pow')
ylabel('rhm 6-12Hz pow')
title('ncp X rhm, JPDF')


figure,
imagesc(linspace(-pi,pi,15),pedges(2,:),apr),axis xy
caxis([0,0.12])
title('Filtered[5-13] RHM peak NCP phase Distrb')
xlabel('Phase (radians)')
ylabel('log10 NCP 5-13hz Power')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
state = 'w';

sind = Trial.stc.states{s};


ind = nniz(vh(sind));
vhs = abs(log10(vh(ind)));
vhlim =prctile(vhs,[5,98]);
edges = linspace(vhlim(1),vhlim(2),9);
edges = [edges(1:end-1);edges(2:end)];
edges_labels = mat2cell(10.^mean(edges),1,ones(1,size(edges,2)))';
yss = ys(sind,:,:,:);
yss = yss(ind,:,:,:);
clear ind;

vsc = [];
for i = edges,
    ind = i(1)>=vhs&vhs<i(2);% & 5>=mean(log10(yss(:,fs>13&fs>5,3,3)),2);
    for j = 1:size(yss,3),
        for k = 1:size(yss,3),
            if sum(ind)~=0,
                if k~=j,
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))./mean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
                else
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(yss(ind,:,k,j));
                end

            else
                vsc(find(i(1)==edges(1,:)),:,:) = zeros([1,size(yss,2),size(yss,3)]);
            end

        end
    end
end


for j = 1:size(yss,3),
    for k = 1:size(yss,3),
        subplot2(size(yss,3),size(yss,3),j,k),
        if j~=k
            imagesc(1:size(edges,2),fs,vsc(:,:,j,k)'),axis xy,
            %title([chan_labels{j} ' & ' chan_labels{k} ' coherence'])
        else
            imagesc(1:size(edges,2),fs,log10(vsc(:,:,j,k))'),axis xy,
            %title([chan_labels{j} ' PSD X Body Speed'])
        end
        set(gca,'XTickLabel',cellfun(@num2str,cellfun(@transpose,edges_labels(2:2:8),'UniformOutput',false),'UniformOutput',false)');
        if j~=k,caxis([.3,.9]),end
        if j==size(yss,3),xlabel('Binned Body Speed cm/s'),end
        if k==1,ylabel('Frequency Hz'),end
        colorbar
    end
end
%suptitle(['Coherence During ' sind.label]);

reportfig(Trial, ...
          'FileName','mean_Coherence_RHM_NCP','Comment',['State: ' ...
                    Trial.stc.states{s}.label]  )

end

tt = 1040;
figure,
h = line(xyz.data(tt,[5,10],1),xyz.data(tt,[5,10],2),xyz.data(tt,[5,10],3)),set(h,'color',[0,1,0])
line(xyz.data(tt,[6,10],1),xyz.data(tt,[6,10],2),xyz.data(tt,[6,10],3))
line(xyz.data(tt,[7,10],1),xyz.data(tt,[7,10],2),xyz.data(tt,[7,10],3))
line(xyz.data(tt,[8,10],1),xyz.data(tt,[8,10],2),xyz.data(tt,[8,10],3))
line(xyz.data(tt,[8,10],1),xyz.data(tt,[8,10],2),xyz.data(tt,[8,10],3))
rhm_mar = 5,12,8,13,14,15;
%h = line([xyz.data(tt,[10],1),vMar(tt,1,1)]',[xyz.data(tt,[10],2),vMar(tt,1,2)]',[xyz.data(tt,[10],3),vMar(tt,1,3)]'),set(h,'color',[1,0,0])
%h = line([xyz.data(tt,[10],1),brx(tt,1,1)]',[xyz.data(tt,[10],2),brx(tt,1,2)]',[xyz.data(tt,[10],3),brx(tt,1,3)]'),set(h,'color',[1,0,0])

for i = 1:5;
h = line(xyz.data(tt,[10,11+i],1),xyz.data(tt,[10,11+i],2),xyz.data(tt,[10,11+i],3)),set(h,'color',[1,0,0])
end






% phase difference vs mean ncp[6,12] power

