function drz = pfDRD(xyz,pfcenters)

npc = size(pfcenters,1);


drz = zeros([xyz.size(1),npc]);
for i = 1:npc,
    %head position
    %pfrs = nan([xyz.size(1),1]);
    %pfrv = nan([xyz.size(1),1]);
    if xyz.size(3)>2,
        pfhxy= cat(2,xyz.data,permute(repmat(pfcenters(i,:),[xyz.size(1),1,1]),[1,3,2]));
    else
        pfhxy = cat(2,cat(3,xyz.data(:,:,[1,2]),zeros([xyz.size([1,2]),1])),...
                      permute(repmat([pfcenters(i,[1,2]),0],xyz.size(1),1),[1,3,2]));
    end
    pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
    
    %head vector

    cor = cell(1,3);
    [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
    cor = cell2mat(cor);
    
    %head to place field center vector
    por = cell(1,3);
    [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    por = cell2mat(por);

    pfhxy.filter(gtwin(1,pfhxy.sampleRate));
    spor = cell(1,3);
    [spor{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    spor = cell2mat(spor);

    drz(:,i) = spor(:,3);
    drz(:,i) = drz(:,i).*[-1;sign(diff(Filter0(gtwin(5,xyz.sampleRate),spor(:,3))))];
end
% $$$     pfds = circ_dist(cor(:,1),por(:,1));
% $$$     pfps = circ_dist(cor(:,2),por(:,2));
% $$$ 
% $$$     tang = acos(dot(sq(diff(pfhxy(:,[1,2],:),1,2)),sq(diff(pfhxy(:,[1,3],:),1,2)),2)./(sqrt(sum(sq(diff(pfhxy(:,[1,2],:),1,2)).^2,2)).*sqrt(sum(sq(diff(pfhxy(:,[1,3],:),1,2)).^2,2))));
% $$$     tangents = cat(3,sin(pfds).*cos(pfps), sin(pfds).*sin(pfps),cos(pfds));
% $$$     mxyz = repmat(pfhxy(:,2,:)-pfhxy(:,1,:),[1,size(tangents,2),1]);
% $$$     tpfds = acos(dot(tangents,mxyz,3)./(sqrt(sum(tangents.^2,3)).*sqrt(sum(mxyz.^2,3))));
% $$$ 
% $$$     
% $$$     drz = nan([xyz.size(1),npc]);
% $$$     pfrs(por(:,3)<thresh,i) = por(por(:,3)<thresh,3);
% $$$     cpc = LocalMinima(pfrs,300);
% $$$ 
% $$$     for c = 1:numel(cpc);
% $$$        drz(find(isnan(pfrs(1:cpc(c))),1,'last'):cpc(c),i) =  -por(find(isnan(pfrs(1:cpc(c))),1,'last'):cpc(c),3);
% $$$        drz(cpc(c):cpc(c)+find(isnan(pfrs(cpc(c):end)),1,'first'),i) =  por(cpc(c):cpc(c)+find(isnan(pfrs(cpc(c):end)),1,'first'),3);
% $$$     end
% $$$ 
% $$$     pfrv = Filter0(gausswin(21)./sum(gausswin(21)),diff(Filter0(gausswin(121)./sum(gausswin(121)),drz(:,i)),1));
% $$$     pfrv(pfrv<2) = nan;
% $$$     tper = Trial.stc{'w'};
% $$$     
% $$$     vpc = LocalMinima(-pfrv,300);
% $$$ 
% $$$     tper.data = [vpc-1,vpc+1];
% $$$     tper = tper+[-2,2];
%end

% $$$ pdst = MTADxyz('data',por(:,3),'sampleRate',xyz.sampleRate);
% $$$ ptang = MTADxyz('data',-sign(tang-pi/2),'sampleRate',xyz.sampleRate);
% $$$ tper = resample(Trial.stc{'t'}.cast('TimeSeries'),xyz);
% $$$ 

% $$$ 
% $$$ figure,
% $$$ sp(1) = subplot(211);
% $$$ plot(pdst(:)),hold on,plot(find(tper.data==1),pdst(tper(:)==1),'.g')
% $$$ sp(2) = subplot(212);
% $$$ plot(ptang(:)),ylim([-3,3])
% $$$ hold on
% $$$ plot(find(tper.data==1),ptang(tper(:)==1),'.g')
% $$$ linkaxes(sp,'x');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ ang = Trial.ang.copy;ang.create(Trial,xyz);
% $$$ qiv = cell(1,3);
% $$$ [qiv{:}] = sph2cart(ang(:,5,7,1),ang(:,5,7,2),ones([ang.size(1),1]).*100);
% $$$ qiv = cell2mat(qiv);
% $$$ 
% $$$ 
% $$$ ind = 118500:10:119000;
% $$$ ind = 417500:10:419000;
% $$$ figure,hold on
% $$$ plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.')
% $$$ scatter3(xyz(ind(1:3),7,1),xyz(ind(1:3),7,2),xyz(ind(1:3),7,3),20,'g')
% $$$ scatter3(xyz(ind(end-3:end),7,1),xyz(ind(end-3:end),7,2),xyz(ind(end-3:end),7,3),20,'r')
% $$$ scatter3(pfcenters(i,1),pfcenters(i,2),pfcenters(i,3),20,'k')
% $$$ quiver3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),qiv(ind,1),qiv(ind,2),qiv(ind,3),0)
% $$$ xlim([-500,500]);
% $$$ ylim([-500,500]);
% $$$ zlim([0,300]);

