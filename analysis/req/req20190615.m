% req20190615
%    Tags: IF04 data processing
%    Status: active
%    Type: analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: NA
%    Description: generate behvaioral states from only head movements for gerrit

% protocol:
%     1. DEFINE states for classification
%     3. DEFINE feature set for each state
%     4. FIND state periods using feature sets
%     5. ENCAPSULATE in MTADepoch object



% LOAD data
Trial = MTATrial.validate('IF04-20181012b.cpm.all');
Trial = MTATrial.validate('IF04-20181006a.spb.all');



% DEFINE states for classification
% state -- description -------------------------------------
% rear  | elevation of upper body over the hindlimbs 
% hdip  | head dip
% loc   | translation through
% pause | immobile while in prone posture (not laying) 
% sit   | immobile while in prone posture (laying down)
stateNames  = {'rearing','head dip','locomotion','turn','pause','sit','rest area','decision point'};
stateLabels = {   'rear',    'hdip',       'loc','turn','pause','sit',     'rest',        'dpoint'};
stateKeys   = {      'r',       'd',         'x',   'n',    'p',  's',        'y',             'i'};




% DEFINE feature set for each state
% state -- feature
% rear  | { head COM height; positive pitch @ onset; positive acceleration @ onset; negative acceleration @ offset }
% hdip  | { head COM height; negative head pitch }
% loc   | { forward head speed high } 
% pause | { forward head speed low }
% sit   | { forward head speed very low }

xyz = Trial.load('xyz','trb');
xyz.filter('ButFilter',3,1.5);

vxy = xyz.vel('head_neck',[1,2]);
vxy.path = xyz.path;
vxy.label = 'headSpeed';
vxy.name = 'head speed';
vxy.ext = 'fet';
vxy.sync = xyz.sync.copy();
vxy.origin = xyz.origin;
vxy.updateFilename(Trial);
vxy.save();

lvxy = vxy.copy();
lvxy.data(vxy.data>0) = log10(vxy.data(vxy.data>0));

rhm = fet_rhm(Trial);
rhm.updateFilename(Trial);
rhm.save();

hrf = fet_href_H(Trial);
fhrf = filter(hrf.copy(),'ButFilter',3,2.5,'low');

[U,mu,vars] = pca(hrf(nniz(xyz),:));

figure,
plot(mu)
hold('on')
plot(rhm.data*10);

uhrf = hrf.copy();
uhrf.data(nniz(xyz),:) = U.*6000;
uhrf.data(nniz(xyz),3) = mu;


fuhrf = filter(uhrf.copy(),'ButFilter',3,1.5,'low');
fuhrf.label = 'fet_href_H_pca_filtLow';
fuhrf.updateFilename(Trial);
fuhrf.save(); 


figure,plot(mu,



figure,
    hold('on');
    plot(vxy(:,1)/2);
    plot(uhrf(:,1:2).*6000);
    plot(fuhrf(:,1:2).*6000);
    Lines([],2.5,'k');
    Lines([],-15,'k');
    Lines([],15,'k');    


lfuhrf = filter(uhrf.copy(),'ButFilter',3,0.5,'low');
figure,
    hold('on');
    plot(vxy(:,1)/2);
    plot(uhrf(:,1:2).*6000);
    plot(lfuhrf(:,1:2).*6000);
    Lines([],2.5,'k');
    Lines([],-15,'k');
    Lines([],15,'k');    
    
figure,
    hold('on');
    plot(xyz(:,'head_neck',1),fuhrf(:,2).*6000,'.');
    Lines([],0,'k');

figure,
    plot(xyz(:,'head_neck',1),xyz(:,'head_neck',2),'.')

    
ang = create(MTADang,Trial,xyz);


figure,
plot(ang(:,'head_back','head_front',2),xyz(:,'head_neck',3),'.');


figure,
    histogram2(ang(:,'head_back','head_front',2),...
               xyz(:,'head_neck',3),...
               [150,150],...
               'XBinLimits',[-pi/2,pi/2],...
               'YBinLimits',[-50,160],...
               'DisplayStyle','tile');

figure,
    histogram2(fuhrf(:,1).*6000,...
               xyz(:,'head_neck',3),...
               [150,150],...
               'XBinLimits',[-10,100],...
               'YBinLimits',[80,200],...
               'DisplayStyle','tile');

figure,
    histogram2(lvxy(:,1),...
               xyz(:,'head_neck',3),...
               [150,150],...
               'XBinLimits',[-3,2],...
               'YBinLimits',[-50,160],...
               'DisplayStyle','tile');

% 'rear','hdip','loc','turn','pause','sit','rest';    
%         {                         fet;      thresh;        fun;       minPer   };
svars = { {        xyz(:,'head_neck',3);         120;        @gt;          30    },  ... % rear
          {        xyz(:,'head_neck',3);           0;        @lt;          30    },  ... % head dip          
          {            fuhrf(:,1).*6000;         2.5;        @gt;          30    },  ... % locomation
          {       abs(fuhrf(:,2).*6000);          15;        @gt;          30    },  ... % turn
           cat(2,{abs(fuhrf(:,1).*6000);         2.5;        @lt;          30    },  ... % pause               
                 {abs(fuhrf(:,2).*6000);          15;        @lt;          30    },  ...
                 { xyz(:,'head_neck',3);          30;        @gt;          30    }), ... 
           cat(2,{ xyz(:,'head_neck',3);           0;        @gt;         300    },  ... % sit
                 { xyz(:,'head_neck',3);          30;        @lt;         300    },  ...
                 {            lvxy(:,1);       -0.75;        @lt;         300    }), ...
           cat(2,{ xyz(:,'head_neck',2);        -450;        @lt;         300    },  ... % rest area
                 {abs(xyz(:,'head_neck',1));     400;        @lt;         300    }), ...
           cat(2,{ xyz(:,'head_neck',2);         450;        @gt;         300    },  ... % decision point
                 {abs(xyz(:,'head_neck',1));     250;        @lt;         300    })  ...           
          
};


Stc = Trial.stc.copy();
Stc.states = {};

for s = 1:numel(svars),
    state = true([size(xyz,1),1]);
    for j = 1:size(svars{s},2),        
        state = state & svars{s}{3,j}(svars{s}{1,j},svars{s}{2,j});
    end

    statePeriods = ThreshCross(state,0.5,svars{s}{4,j});

    Stc.addState(Trial.spath,...
                 Trial.filebase,...
                 statePeriods,...
                 xyz.sampleRate,...
                 Trial.sync.copy,...
                 Trial.sync.data(1),...
                 stateLabels{s},stateKeys{s});
end
    
Trial.stc = Stc;
Trial.stc.save(1);
Trial.save;