BodyRef = [];
for ii=2:5;
    
   tmp= bsxfun (@minus,Trial.xyz.data(:,ii,:),Trial.xyz.data(:,ii-1,:));
   Ref = bsxfun(@rdivide,tmp,sqrt(sum(tmp.^2,3))); %unit vector notation
    distance=[];
distance=(sqrt(nansum(tmp(:,1,:).^2,3)));
distance(distance==0)=nan;
distanceNorm= (distance-nanmean(distance))/nanvar(distance);

[Az,Pit,xx]=cart2sph(Ref(:,1,1),Ref(:,1,2),Ref(:,1,3));

BodyRef(:,:,ii-1)=[Az, Pit ,distanceNorm];

end
tmp=[];

AngDiffBody=[];
for ii =2:4;

AngDiffBody(:,1,ii-1)= circ_dist(BodyRef(:,1,ii-1),BodyRef(:,1,ii));
AngDiffBody(:,2,ii-1)= circ_dist(BodyRef(:,2,ii-1),BodyRef(:,2,ii));
AngDiffBody(:,3,ii-1)= nanvar([BodyRef(:,3,ii-1),BodyRef(:,3,ii)],0,2);
end


% Az1=(circ_dist(AngDiffBody(:,1,1),AngDiffBody(:,1,2)));
% 
% Az2=(circ_dist(AngDiffBody(:,1,2),AngDiffBody(:,1,3)));
% 
% Pit1=(circ_dist(AngDiffBody(:,2,1),AngDiffBody(:,2,2)));
% 
% Pit2=(circ_dist(AngDiffBody(:,2,2),AngDiffBody(:,2,3)));
% 
% 
% AzDiff = circ_dist(Az1,Az2);
% AzDiff(isnan(AzDiff))=0;
% PiDiff = circ_dist(Pit1,Pit2);
% PiDiff(isnan(PiDiff))=0;


%%
%distance
DistVar=[];
tmp =sum(AngDiffBody(:,3,:),3);
% DistVar = (tmp-nanmean(tmp))/nanvar(tmp);
DistVar = 10*((tmp-min(tmp))/max(tmp)-min(tmp));

%Az
X = sin(sq(AngDiffBody(:,1,:)));
Y = cos(sq(AngDiffBody(:,1,:)));

Xm = nanmean(X,2);
Ym = nanmean(Y,2);
R=[];
R= sqrt(sum([Xm.^2 Ym.^2],2));
% R = (R-nanmean(R))/nanvar(R);
RV = R-DistVar;
mAng = atan2(Ym,Xm);

%Pit

Xp = sin(sq(AngDiffBody(:,2,:)));
Yp = cos(sq(AngDiffBody(:,2,:)));

Xmp = nanmean(Xp,2);
Ymp = nanmean(Yp,2);

Rp=[];
Rp= sqrt(nansum([Xmp.^2 Ymp.^2],2));
% Rp = (Rp-nanmean(Rp))/var(Rp);
RVp = Rp-DistVar;
mAngp = atan2(Ymp,Xmp);

medR=medfilt1(RV,round(1*Trial.xyz.sampleRate));
medRp=medfilt1(RVp,round(1*Trial.xyz.sampleRate));


TestErrors = sum([medR<.8 medRp<.7 Trial.xyz.data(:,1,3)<0],2)>0;


figure
imagesc([],-2:2,[Errors1 TestErrors*2]'); axis xy; colormap 'parula'
hold on
plot(.005*(Trial.xyz.data(:,7,3)),'m','LineWidth',2)
plot(.005*(Trial.xyz.data(:,1,3)),'g','LineWidth',2)

%%
[Errors1,Errors4m]=ExtractErrors(Trial,2);

figure
imagesc([],-5:5,[~Errors1]');hold on;axis xy; colormap 'gray'


Test= ButFilter((circ_dist(Az1,Az2)),[4],[1/(Trial.xyz.sampleRate/2) ],'low' );


TestFilt(isnan(TestFilt))=0;

TestFilt1= ButFilter(unwrap(TestFilt),[4],[.01/(Trial.xyz.sampleRate/2) ],'high');


figure
imagesc([],-3:3,(~(TestFilt<0.4))'); axis xy; colormap 'gray'
hold on
plot((circ_dist(Az1,Az2)))
plot((circ_dist(Pit1,Pit2)))
plot(TestFilt,'r')


%%
figure

imagesc([],-5:5,[((TestFilt<0.4)) Errors1*2]'); axis xy; 
hold on
plot(.03 * (Trial.xyz.data(:,1,3)),'g')
plot(.03 * (Trial.xyz.data(:,7,3)),'m')




