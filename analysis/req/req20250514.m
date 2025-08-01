
% req20250514
% Using the spike features and accg to sort units.


bpc_load_data();

cf(@(t) t.load('nq'), Trials);
ampSymInt = cf(@(T,U) T.nq.AmpSym(U), Trials,UnitsInt);
ampSymInt = cat(1,ampSymInt{:});
ampSymPyr = cf(@(T,U) T.nq.AmpSym(U), Trials,Units);
ampSymPyr = cat(1,ampSymPyr{:});

spkWidthRInt = cf(@(T,U) T.nq.SpkWidthR(U), Trials,UnitsInt);
spkWidthRInt = cat(1,spkWidthRInt{:});
spkWidthRPyr = cf(@(T,U) T.nq.SpkWidthR(U), Trials,Units);
spkWidthRPyr = cat(1,spkWidthRPyr{:});

centerMaxInt = cf(@(T,U) T.nq.CenterMax(U), Trials,UnitsInt);
centerMaxInt = cat(1,centerMaxInt{:});
centerMaxPyr = cf(@(T,U) T.nq.CenterMax(U), Trials,Units);
centerMaxPyr = cat(1,centerMaxPyr{:});

Accg = cf(@(T) autoccg(T), Trials);
Accg = cat( 2, Accg{:});
Accg_O = Accg;
Accg = Accg_O;

ampSym = cf(@(T) T.nq.AmpSym, Trials);
ampSym = cat(1,ampSym{:});

spkWidthR = cf(@(T) T.nq.SpkWidthR, Trials);
spkWidthR = cat(1,spkWidthR{:});

figure();
imagesc(log10(Accg)')
figure();
histogram(log10(sum(Accg)));

spkWidthR(sum(Accg)<400)=[];
ampSym(sum(Accg)<400)=[];
Accg(:,sum(Accg)<400)=[];




Nccg = Accg;
Nccg = bsxfun(@rdivide,Nccg,sum(Nccg));


[LU, LR, FSr, VT] = erpPCA( Nccg', 4, 100 );

figure(); plot(VT(1:4,1),'-+');

figure();
imagesc(LU');

figure();
plot(FSr(:,1),FSr(:,3),'.');
figure();
plot(FSr(:,2),FSr(:,3),'.');

FSrO = FSr;

FSr(FSr(:,1)<-5,:) = [];



FSrXY = tsne([FSr,ampSym,spkWidthR],[],2,5:6,40);


figure();
plot(FSrXY(:,1),FSrXY(:,2),'.');

figure();
subplot(3,3,1);
plot(FSrXY(:,1),FSrXY(:,2),'.');
g1 = ~(FSrXY(:,2) > 0 & (FSrXY(:,1) < 900 | FSrXY(:,1) > 40));
g2 = FSrXY(:,2) > 0 & (FSrXY(:,1) < 900 | FSrXY(:,1) > 40);
subplot(3,3,2); imagesc(log10(Nccg(:, g1)'));
subplot(3,3,3); imagesc(log10(Nccg(:, g2)'));
subplot(3,3,4); imagesc(log10(Nccg(:, FSrXY(:,1) < -80 & (FSrXY(:,1) < 90 | FSrXY(:,1) > -90))'));
subplot(3,3,5); imagesc(log10(Nccg(:, FSrXY(:,2) < -30 & (FSrXY(:,1) < -20 | FSrXY(:,1) > -60))'));
subplot(3,3,6); imagesc(log10(Nccg(:, FSrXY(:,2) > 20 & (FSrXY(:,1) < 90 | FSrXY(:,1) > -20))'));
subplot(3,3,7); imagesc(log10(Nccg(:, FSrXY(:,2) > -20 & (FSrXY(:,1) < -20 | FSrXY(:,1) > -70))'));



[LU1, LR1, FSr1, VT1] = erpPCA( Nccg(:,g1)', 3, 100 );
figure,plot(VT1(:,1:3))
figure,plot(FSr1(:,1),FSr1(:,2),'.');

FSrXY1 = tsne([FSr1],[],2,1:2,100);
figure();
plot(FSrXY1(:,1),FSrXY1(:,2),'.');

[~,sind] = sort(FSrXY1(:,1));

Nccg1 = Nccg(:,g1);
figure();
imagesc(log10(Nccg1(:,sind)'));


[LU2, LR2, FSr2, VT2] = erpPCA( Accg(:,g2)', 4, 100 );
figure,plot(VT2(:,1))
figure,plot(FSr2(:,1),FSr2(:,2),'.');



g2(g2) = FSr2(:,2)<10;

FSrXY2 = tsne([FSr2],[],2,2:4,20);
figure();
plot(FSrXY2(:,1),FSrXY2(:,2),'.');
g21 = g2;
g21(g21) = FSrXY2(:,1)>50;
g22 = g2;
g22(g22) = FSrXY2(:,1)<50;
g23 = g2;
g23(g23) = FSrXY2(:,1)<15 & FSrXY2(:,2)>40;

[~,sind] = sort(FSrXY2(:,1));
figure();
imagesc(log10(Nccg(:,g21)'));
figure();
imagesc(log10(Nccg(:,g22)'));
figure();
imagesc(log10(Nccg(:,g23)'));


Nccg2 = Nccg(:,g2);
[~,sind] = sort(Nccg2(53,:));
figure();
imagesc(log10(Nccg2(:,sind)'));

[~,sind] = sort(Nccg1(69,:));
[~,sind] = sortrows(Nccg1([46,69],:)', 'ascend');
figure();
imagesc(log10(Nccg1(:,sind)'));



figure();
subplot(3,3,1);
plot(FSrXY(:,1),FSrXY(:,2),'.');
subplot(3,3,2); imagesc(log10(Nccg(:, FSrXY(:,1) < 60 & (FSrXY(:,1) < 0 | FSrXY(:,1) > -105))'));
subplot(3,3,3); imagesc(log10(Nccg(:, FSrXY(:,2) < 0 & (FSrXY(:,1) < 60 | FSrXY(:,1) > -20))'));
subplot(3,3,4); imagesc(log10(Nccg(:, FSrXY(:,1) < -80 & (FSrXY(:,1) < 90 | FSrXY(:,1) > -90))'));
subplot(3,3,5); imagesc(log10(Nccg(:, FSrXY(:,2) < -30 & (FSrXY(:,1) < -20 | FSrXY(:,1) > -60))'));
subplot(3,3,6); imagesc(log10(Nccg(:, FSrXY(:,2) > 20 & (FSrXY(:,1) < 90 | FSrXY(:,1) > -20))'));
subplot(3,3,7); imagesc(log10(Nccg(:, FSrXY(:,2) > -20 & (FSrXY(:,1) < -20 | FSrXY(:,1) > -70))'));
