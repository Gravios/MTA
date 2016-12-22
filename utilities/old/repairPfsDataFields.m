%function dstruct = population_placeFieldAnalysis(sesList,states,pftype)

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','i_rt','i_wt','i_gt','i_lt'};

pftype = 'MTAAknnpf';

numsts = numel(states);
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});

    for i = 1:numsts,
        switch pftype
          case 'MTAAknnpf'
            filebase = {sesList{ses}{1},'.',sesList{ses}{2},'.',sesList{ses}{3}, ...
                        '.pfknn.xy.head_front.',states{i}, '.us0.5bs1000nnn80dt70bd20_20ds20.mat'};
            filebase = strjoin(filebase,'');

            if exist(fullfile(Trial.spath,filebase),'file')
                load(fullfile(Trial.spath,filebase))
            else
                continue,
            end
            
            field2rm = {'maxRate',...
                        'maxRateInd',...
                        'maxRatePos',...
                        'meanRate',...
                        'si',...
                        'spar'};
            try
            for j = 1:numel(field2rm),
                Pfs.data = rmfield(Pfs.data,field2rm{j});
            end
            end

            Pfs.path = Trial.spath;
            save(Pfs.fpath,'Pfs');
          case 'MTAApfs'
            % not ready 
            % pfs{i} = MTAApfs(Trial,units,states{i},0,'numIter',1000)
        end

    end
end