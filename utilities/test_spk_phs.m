Trial = MTATrial('jg05-20120315');
cd(Trial.spath);
[Thph,ThAmp,ThFr] = ThetaParams(Trial.name,61,1,[],[],[],[]);

[Res,Clu,Map] = LoadCluRes(Trial.name);
Res = round(Res./Trial.sampleRate.*Trial.lfp.sampleRate);



tper = load(fullfile(Trial.spath,[Trial.name '.sts.theta']));

myRes = Res(Clu==32);
sResTheta = SelectPeriods(myRes,tper,'d',1,1);
ThphTheta = SelectPeriods(Thph,tper,'c',1);

figure,
circ_plot(ThphTheta(sResTheta),'hist',[],30,true,true);

