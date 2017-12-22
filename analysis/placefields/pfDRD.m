function drz = pfDRD(xyz,pfs)


[~,pfCenters] = pfs.maxRate;

npc = size(pfCenters,1);

drz = zeros([xyz.size(1),npc]);

for i = 1:npc,
    if xyz.size(3)>2,
        pfhxy= cat(2,xyz.data,permute(repmat(pfCenters(i,:),[xyz.size(1),1,1]),[1,3,2]));
    else
        pfhxy = cat(2,cat(3,xyz.data(:,:,[1,2]),zeros([xyz.size([1,2]),1])),...
                      permute(repmat([pfCenters(i,[1,2]),0],xyz.size(1),1),[1,3,2]));
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
