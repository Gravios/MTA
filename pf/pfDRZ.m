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
    else
        pfrs(:,i) = por(:,3).*cos(por(:,2));
    end
    pfds(:,i) = circ_dist(cor(:,1),por(:,1));
    pfps(:,i) = circ_dist(cor(:,2),por(:,2));
    
end

if xyz.size(3)>2,
tangents = cat(3,sin(pfds).*cos(pfps), sin(pfds).*sin(pfps),cos(pfds));
mxyz = repmat(pfhxy(:,2,:)-pfhxy(:,1,:),[1,size(tangents,2),1]);
tpfds = acos(dot(tangents,mxyz,3)./(sqrt(sum(tangents.^2,3)).*sqrt(sum(mxyz.^2,3))));
pfds = zeros(size(tpfds));
pfds(abs(tpfds)<=pi/2)=1;
pfds(abs(tpfds)>pi/2)=-1;
else
    tpfds = pfds;
    pfds(abs(tpfds)<=pi/2)=-1;
    pfds(abs(tpfds)>pi/2)=1;
end

%DRZ
drz = pfds.*pfrs;
