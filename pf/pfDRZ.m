function drz = pfDRZ(Trial,xyz,pfs,varargin)

[dims] = DefaultArgs(varargin,{[1,2]},1);

[pfMaxRates,pfCenters] = pfs.maxRate;
pfMaxRates = repmat(pfMaxRates(:)',xyz.size(1),1);

numPfs = size(pfCenters,1);

pfrs = zeros([xyz.size(1),numPfs]);
pfps = zeros([xyz.size(1),numPfs]);
pfds = zeros([xyz.size(1),numPfs]);


for i = 1:numPfs,
    %head position
    pfhxy = xyz.copy;
    if dims>2,
        pfhxy.data= cat(2,xyz.data,permute(repmat(pfCenters(i,:),[xyz.size(1),1,1]),[1,3,2]));
    else
        pfhxy.data = cat(2,cat(3,xyz.data(:,:,dims),zeros([xyz.size([1,2]),1])),...
                      permute(repmat([pfCenters(i,dims),0],xyz.size(1),1),[1,3,2]));
    end
    
    %head vector
    cor = cell(1,3);
    [cor{:}] = cart2sph(pfhxy(:,end,1)-pfhxy(:,'',1),...
                        pfhxy(:,end,2)-pfhxy(:,1,2),...
                        pfhxy(:,end,3)-pfhxy(:,1,3));
    cor = cell2mat(cor);
    
    %head to place field center vector
    por = cell(1,3);
    [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    por = cell2mat(por);
    
    % Calculate head-placefield orientation (pfds)
    % and head-placefield distance (pfrs)
    if xyz.size(3)>2,
        pfds(:,i) = sign(acos(dot(sq(diff(pfhxy(:,[1,2],:),1,2)),sq(diff(pfhxy(:,[1,3],:),1,2)),2)./(sqrt(sum(sq(diff(pfhxy(:,[1,2],:),1,2)).^2,2)).*sqrt(sum(sq(diff(pfhxy(:,[1,3],:),1,2)).^2,2))))-pi/2);    
    else
        pfds(:,i) = sign(circ_dist(cor(:,1),por(:,1))-pi/2);
    end

end


%DRZ
drz = pfds.*pfds;
