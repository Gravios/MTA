function fet = detect_walk(Trial)

xyz = Trial.xyz.copy; 
xyz.load(Trial);
xyz.filter(gausswin(91)./sum(gausswin(91)));

fet = [];
w = round([0.5,1,1.5,2]).*xyz.sampleRate;
for m = 2:5,
    for k = 1:numel(w),
        wins = round(linspace(10,w(k),30));
        nw = numel(wins);
        sxy = zeros([xyz.size([1,2]),2,nw]);
        for i = 1:9,
            for j = 1:nw
                sxy(:,i,:,j) = circshift(xyz(:,i,[1,2])-circshift(xyz(:,i,[1,2]),-wins(j)),round(wins(j)/2));
            end
        end
        ims = circshift(sq(sum(sxy(:,3,:,1:30).*sxy(:,m,:,1:30),3)),round(-.5*wins(end)))';
        fet(:,k,m-1) = mean(diff(diff(ims')'))';
    end
end

fet = sq(mean(fet,2));
fet = sq(mean(fet,2));
fet = [0;fet];
