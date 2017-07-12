% rf_ccg_grps.m
% No fucking clue
% maybe something to do with catching patch sizes?

load ~/data/analysis/jg05-20120317/jg05-20120317.cof.all.ccg.rear_pfc_walk.mat

frccg = Bccg.filter(gausswin(3));


figure,
subplot2(2,1,1,1); hist(sq(max(frccg(:,51,3,1,:),[],1)),200),
subplot2(2,1,2,1); hist(sq(min(frccg(:,51,3,1,:),[],1)),200),


sccg = sort(sq(frccg(:,51,3,1,:)),2);

sthr = sccg(:,round(size(sccg,2)*.95))/2;



for i = 1:20,
figure,
hold on,
Bccg.tbin,frccg(:,51,3,1,i));
plot(Bccg.tbin',sthr,'r');
end


clear('sccg')
clear('sthr')

sccg = sort(frccg,5);
sthr = repmat(sccg(:,:,3,:,round(size(sccg,5)*.95))/2,[1,1,3,1,size(sccg,5)]);
sthr = repmat(max(sccg(:,:,3,:,round(size(sccg,5)*.95)))/2,[size(sccg,1),1,3,1,size(sccg,5)]);

grps = false(size(frccg));
grps(sccg>sthr)=true;

ccgsz=size(sccg);
gsta = zeros([round(size(sccg,1)/2*size(sccg,5)),ccgsz(2:4),2]);

b=3;g=1;k=1;
for u = 1:size(sccg,2),
    for p = 1:4,
g = 1;
        for i = 1:size(sccg,5)
            k=1;
            gind = find(grps(:,u,b,p,i));
            if isempty(gind),continue,end
            gnel = numel(gind);
            gend=0;
            while k<gnel
                gdif = diff(gind(k:end));
                gend = find(gdif~=1,1,'first');
                if isempty(gend)&~isempty(gdif),
                    gend = numel(gdif)+1;
                elseif isempty(gend)&isempty(gdif),
                    gend = gnel+1;
                end
                gsta(g,u,b,p,1) = max(sccg(gind(k):gind(k)+gend-1,u,b,p,i));
                gsta(g,u,b,p,2) = gend;
                g=g+1;
                if k>gnel,
                    k=1;
                    continue,
                else
                    k=find(gind==(gind(k)+gend-1),1,'first')+1;
                end
            end
        end
    end
end

b=2;g=1;k=1;
for u = 1:size(sccg,2),
    for p = 1:4,
g = 1;

            k=1;
            gind = find(grps(:,u,b,p,1));
            if isempty(gind),continue,end
            gnel = numel(gind);
            gend=0;
            while k<gnel
                gdif = diff(gind(k:end));
                gend = find(gdif~=1,1,'first');
                if isempty(gend)&~isempty(gdif),
                    gend = numel(gdif)+1;
                elseif isempty(gend)&isempty(gdif),
                    gend = gnel+1;
                end
                gsta(g,u,b,p,1) = max(sccg(gind(k):gind(k)+gend-1,u,b,p,1));
                gsta(g,u,b,p,2) = gend;
                g=g+1;
                if k>gnel,
                    k=1;
                    continue,
                else
                    k=find(gind==(gind(k)+gend-1),1,'first')+1;
                end
            end

    end
end




[acccg,tbin] = autoccg(Session);
pfr = MTAPlaceField(Session,[],'rear');
pfw = MTAPlaceField(Session,[],'walk');

unit =1;
while unit~=-1,
clf
subplot2(3,5,1,1);
pfr.plot(unit,1);
subplot2(3,5,2,1);
pfw.plot(unit,1);
subplot2(3,5,3,1);
bar(tbin,accg(:,unit));axis tight
for p = 2:5,
subplot2(3,5,1,p);,try hist2(snq(gsta(gsta(:,unit,3,p-1,1)~=0,unit,3,p-1,[2,1])),20,20),end
                   hold on,
                   plot(gsta(gsta(:,unit,1,p-1,1)~=0,unit,1,p-1,2),gsta(gsta(:,unit,1,p-1,1)~=0,unit,1,p-1,1),'*w')
                   plot(gsta(gsta(:,unit,2,p-1,1)~=0,unit,2,p-1,2),gsta(gsta(:,unit,2,p-1,1)~=0,unit,2,p-1,1),'*r')
                   cac = caxis;
                   caxis(cac./4)
                   fprintf('unit:%i | caxis:%d %d\n',[unit,cac])
subplot2(3,5,2,p); bar(Bccg.tbin,sccg(:,unit,1,p-1,1));
                   axis tight
                   Lines([],sthr(1,unit,1,p-1,1));
subplot2(3,5,3,p); bar(Bccg.tbin,sccg(:,unit,2,p-1,1));
                   axis tight
                   Lines([],sthr(1,unit,1,p-1,1));
end
title(num2str(unit));
unit = figure_controls(gcf,unit);
end


figure
unit=14
while unit~=-1,
clf
subplot(211)
pfr.plot(unit,1);
subplot(212)
pfw.plot(unit,1);
title(num2str(unit));
unit = figure_controls(gcf,unit);
end

 x=0:10;
 y=x.^2;
 err=0.1*y;
 figure; hold on; plot(x,y,'r.-');
 patch([x fliplr(x)],[y+err fliplr(y-err)],[0.7 0.7 0.7]);
plot(x,y,'r');