function drz = pfDRZ(xyz,pfcenters)

npc = size(pfcenters,1);

pfrs = zeros([xyz.size(1),npc]);
pfps = zeros([xyz.size(1),npc]);
pfds = zeros([xyz.size(1),npc]);
for i = 1:npc,
    %head position
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
    if xyz.size(3)>2,
        pfrs(:,i) = por(:,3);
        pfds(:,i) = sign(acos(dot(sq(diff(pfhxy(:,[1,2],:),1,2)),sq(diff(pfhxy(:,[1,3],:),1,2)),2)./(sqrt(sum(sq(diff(pfhxy(:,[1,2],:),1,2)).^2,2)).*sqrt(sum(sq(diff(pfhxy(:,[1,3],:),1,2)).^2,2))))-pi/2);    
    else
        pfds(:,i) = sign(circ_dist(cor(:,1),por(:,1))-pi/2);
        pfrs(:,i) = por(:,3).*cos(por(:,2));
    end
end


%DRZ
drz = pfds.*pfrs;
