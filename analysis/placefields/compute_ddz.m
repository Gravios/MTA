function ddz = compute_ddz(Trial,pft,units)


pfstats = compute_pfstats_bs(Trial);
dims = 2;

[mrt,mrp] = pft.maxRate(units);
xyz = Trial.load('xyz');

% Get the rat's heading 
pfds = [];
pfdd = [];
for unit = units
    pfhxy = xyz(:,{'head_back','head_front'},:);
    %pfhxy = cat(2,pfhxy,permute(repmat([mrp(unit==units,:),0],[size(xyz,1),1]),[1,3,2]));
    pfhxy = cat(2,pfhxy,permute(repmat([fliplr(sq(mean(pfstats.peakPatchCOM(8,:,pfstats.cluMap==unit,:)))'),0],[size(xyz,1),1]),[1,3,2]));    

    pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);

    cor = cell(1,3);
    [cor{:}] = cart2pol(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
    cor = cell2mat(cor);
    por = cell(1,3);
    [por{:}] = cart2pol(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    por = cell2mat(por);

    
    pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
    pfdd(:,unit==units) = por(:,2);
    
end
pfd = zeros(size(pfds));
pfd(abs(pfds(:))<=pi/2)=-1;
pfd(abs(pfds(:))>pi/2)=1;

% Calculate DRZ 
ddz = pfd.*pfdd;


