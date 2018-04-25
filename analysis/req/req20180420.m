

Trial = MTATrial.validate('jg05-20120311.cof.all');

Trial.sync.data



s = MTASession.validate(Trial.filebase);

xyz = preproc_xyz(Trial,'trb');
ufr = create(MTADufr,s,xyz,[],[],120,true);
ufr.data = bsxfun(@rdivide,ufr.data,prctile(ufr.data,90));

per = [   0,2500;...
       2501,3100;...
       3101, 12669];


shiftPeriods = round(([ufr.origin,2500;...
                       2501,3100;...
                       3100, 12669.360425397]...
                     -ufr.origin).*xyz.sampleRate)+1;

rates = round([sum(ufr(shiftPeriods(1,:),:))'./diff(shiftPeriods(1,:))...
               sum(ufr(shiftPeriods(2,:),:))'./diff(shiftPeriods(2,:)),...
               sum(ufr(shiftPeriods(3,:),:))'./diff(shiftPeriods(3,:))],3);

figure();
bar([0:0.05:1],histc(rates(:,2),[0:0.05:1]),'histc');

thresh = 0.2;

s.spk.per = MTADepoch([],[],          ... Sync Periods
                      shiftPeriods,            ... Data
                      xyz.sampleRate,              ... Sample Rate
                      xyz.sync.sync.sync.copy(),    ... Sync Periods
                      0);             %   Sync Sync Origin

s.spk.perInd = rates>0.2;

s.save

pft = pfs_2d_theta(Trial,[],false,true,1);

figure();
for unit = 1:210,
    plot(pft,unit,1,true,[],true);
    title(num2str(unit));
    pause(0.6);
end

