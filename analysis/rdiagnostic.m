
pf_s = MTAPlaceField([]);
pf_s.stateLabel = 'head.theta';

Session = MTASession('jg05-20120310');
Session = MTATrial(Session,'all');
Session = Session.load_CluRes(Session.xyzSampleRate);
Session = Session.load_Pfs(pf_s);


stsp = round(Session.statePeriods({'theta','head'})./Session.lfpSampleRate.*Session.xyzSampleRate+1);;

myres = SelectPeriods(Session.res(Session.clu==15),stsp,'d',1,1);
mypos = SelectPeriods(Session.xyz(:,7,[1,2]),stsp,'c',1);

figure,plot(mypos(myres,1),mypos(myres,2),'.')

xtsp = round(stsp./Session.lfpSampleRate.*Session.xyzSampleRate+1);
figure,plot(Session.xyz(xtsp,7,1),Session.xyz(xtsp,7,2),'.')

rper = Session.Bhv.getState('rear').state;

%figure,plot(Session.xyz(rper(:,1),7,1),Session.xyz(rper(:,1),7,2),'.')
%hold on,plot(Session.xyz(rper(:,2),7,1),Session.xyz(rper(:,1),7,2),'.','color','r')






grper =  rper(~sum(sum(reshape([Session.xyz(rper,7,1),Session.xyz(rper,7,2)],size(rper,1),2,2)==zeros(size(rper,1),2,2),3),2)>0,:);

[ou,xb,yb,po] = hist2(Session.xyz(Session.xyz(myres,1,1)~=0),7,1:2),50,50);
figure,imagesc(xou'),axis xy

[ous,xbs,ybs,pos] = hist2(Session.xyz(Session.xyz(:,1,1)~=0&sq(sqrt(sum(Session.xyz(:,1,1:2).^2,3)))>600,7,1:2),50,50);

figure,imagesc(xbs,ybs,ous'),axis xy




%% Load Place Fields
Session.Pfs = [];
Session = Session.load_Pfs;


%% Display Everything
figure(102)
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
    clf

    for pf = 1:length(Session.Pfs),
        %% Rate Map
        subplotfit(pf,length(Session.Pfs))
        ppf(Session.Pfs{pf}.xbin,Session.Pfs{pf}.ybin,Session.Pfs{pf}.rateMap{unit})
        title([Session.Pfs{pf}.stateLabel ' ' num2str(unit)])
    end

    %% ReportFig
    %reportfig(102,[Session.filebase '.pfs'],[],['unit: ' num2str(unit)],[],0);
    
    %% Figure controls
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #: ');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
    end
end
