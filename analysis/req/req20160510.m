%function req20160510(Trial)
% $$$ Trial = MTATrial.validate(Trial);
% $$$ 
% $$$ xyz = Trial.load('xyz');
% $$$ fxyz = xyz.filter('ButFilter',3,2.4,'low');
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ vxy = fxyz.vel([],[1,2]);
% $$$ vxy.data(vxy.data<1e-4) = 1e-4;
% $$$ vxy.data = log10(vxy.data);
% $$$ 
% $$$ dang = MTADang('data',circ_dist(ang([end,1:end-1],'pelvis_root','spine_upper',1),...
% $$$                                 ang([1:end],'head_back','head_front',1)),...
% $$$                'sampleRate',xyz.sampleRate);
% $$$ 
% $$$ eds = {linspace(-pi/2,pi/2,100),linspace(-3,2,100)};
% $$$ ind = Trial.stc{'a'};
% $$$ ind = sqrt(sum(xyz(:,7,[1,2]).^2,3))<200;
% $$$ hfig = figure(str2num([cell2mat(regexp('req20160510','\d*','match')),'1']));
% $$$ hist2([dang(ind),vxy(ind,'spine_upper')],eds{:});
% $$$ xlabel('head to body angle (radians)')
% $$$ ylabel('upper spine speed (log10(cm/s)))')
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase],                    ... Tumbnail caption
% $$$           [Trial.filebase,': head body angle within 20cm from center ' ...
% $$$            'vs upper spine speed'],                    ... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$ 

procOpts = {{'newSampleRate',12,'procOpts',''},...
            {'newSampleRate',12,'procOpts','spline_spine'},...
            {'newSampleRate',12,'procOpts','SPLINE_SPINE_HEAD_EQD'}
};
postfix = {'','ss','sshe'};
mapped = {'','-m'};
norm   = {'','-n'};
states = {'rear','walk','turn','pause','groom','sit'};

for s = 1:numel(states),
    for p = 1:3,
        for m = 0:1
            for n = 1
                distributions_features([],...
                                       'fet_mis',...
                                       [],...
                                       ['fet_mis-',postfix{p},mapped{m+1},norm{n+1}],...
                                       'featureOpts',procOpts{p},...
                                       'mapToReference',m,...
                                       'normalize',n,...
                                       'state',states{s})
            end
        end
    end
end
