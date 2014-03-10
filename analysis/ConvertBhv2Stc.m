MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

sname = 'jg05-20120309';
tname = {'all'};
bname = {'m20130311'};
shift = 0;
HL_bhvs = {'rear','walk','shake','turnL','turnR','sniff','groom'};


sname = 'jg05-20120310'; 
tname = {'all'};
bname = {'m20130318'};
shift =  40;
HL_bhvs = {'rear','walk','shake','turnL','turnR','sniff','groom'};


sname = 'jg05-20120315';
tname = {      'crt1',     'alt1',     'crt2',     'alt2',     'crt3'};
bname = {'m20130304','m20130305','m20130307','m20130307','m20130307'};
shift = [       -100,          0,       1000,       1000,       -300];
HL_bhvs = {'rear','walk','shake','turnL','turnR','sniff','groom'};


Trial = MTATrial(sname,'all');
Trial.load('xyz');


dxind = diff(Trial.xyz.sync.data,1,2)
dxind = dxind(1:end-1)+Trial.xyz.sync.data(2:end,1)-Trial.xyz.sync.data(1:end-1,2);
dxind = cumsum(dxind);
dxind = [0;dxind];


Stc = MTAStateCollection(Trial.spath,Trial.filebase,'hand_labeled',Trial.stc.sync.copy,Trial.stc.origin,1,'stc');



for tind = 1:numel(tname),
    load(['/gpfs01/sirota/bach/data/gravio/analysis/' sname '/' sname '.cof.' tname{tind} '.bhv.' bname{tind} '.mat' ]);


    for sind = 1:numel(HL_bhvs),
        bhvind = sind;
        if ~strcmp(Bhv.States{bhvind}.label,HL_bhvs{sind})
            for temp_sind = 1:numel(HL_bhvs),
                if strcmp(Bhv.States{bhvind}.label,HL_bhvs{temp_sind})
                    bhvind = temp_sind;
                    break
                end
            end
        end

        if tind==1,
            Stc.states{sind} = MTADepoch(Trial.spath,...
                                         Trial.filebase,...
                                         (dxind(tind)+shift(tind))+Bhv.States{bhvind}.state,...
                Trial.xyz.sampleRate,...
                Trial.stc.states{1}.sync.copy,...
                Trial.stc.states{1}.origin,...
                Bhv.States{bhvind}.label,...
                Bhv.States{bhvind}.key,...
                'TimePeriods',...
                'sst');
        else
            Stc.states{sind}.data = cat(1,Stc.states{sind}.data,(dxind(tind)+shift(tind))+Bhv.States{bhvind}.state);
        end

    end
end

figure(283848),clf
plot(Trial.xyz(:,7,3))
% LB Rearing Lines
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Stc{'r'}(:),[],'k');
% HL Rearing Lines
%Lines((dxind(tind)+shift(tind))+Bhv.States{2}.state(:),[],'k');






