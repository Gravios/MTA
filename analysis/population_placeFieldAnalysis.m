%function dstruct = population_placeFieldAnalysis(sesList,states,pftype)

sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};

pftype = 'MTAAknnpf';

for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.xyz.load(Trial);
    Trial.load('nq');

    units = [];

    % Load PlaceFields
    pfs={};
    for i = 1:numsts,
        switch pftype
          case 'MTAAknnpf'
            pfs{i} = MTAAknnpfs(Trial,units,states{i},0,'numIter',1000, ...
                            'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
          case 'MTAApfs'
            % not ready 
            % pfs{i} = MTAApfs(Trial,units,states{i},0,'numIter',1000)
        end

        if i==1,
            units = pfs{1}.data.clu;
        end
    end


    % Calulate the characteristic values of place fields
    for u = units
        for state = 1:numsts
            pfstats(state,u) = PlaceFieldStats(Trial,pfs{state},u);    
        end
    end

    % Initialize correlation vars
    stsCor = zeros(numel(units),numsts,numsts,2,2,2);
    % Initialize overlap vars
    overlap = zeros(numel(units),numsts,numsts);
    % Initialize Permutation test vars
    numIter=1000;
    permsamp = zeros(numel(units),numsts,numsts,numIter);

    for u = units,
        for statei = 1:numsts,
            for statej = 1:numsts,

                hsig = pfs{statei}.plot(u,'sig');
                lsig = pfs{statej}.plot(u,'sig');
                hmap = pfs{statei}.plot(u);
                lmap = pfs{statej}.plot(u);
                mind = hsig<.05&lsig<.05;

                if sum(mind(:))>10, 
                    % state ratemap correlations
                    [stsCor(u==units,statei,statej,:,:,1),...
                     stsCor(u==units,statei,statej,:,:,2)] = corrcoef(hmap(mind),lmap(mind));
                    % States ratemap overlap
                    novlp = sum(mind(:))/(sum(reshape(hsig<.05,[],1))+ ...
                                          sum(reshape(hsig<.05,[],1)));
                    overlap(u==units,statei,statej) = novlp;
                    % States permuted mean rate
                    perminds = randi(sum(mind(:))*2,[sum(mind(:)),numIter]);
                    permpool = [hmap(mind);lmap(mind)];
                    permsamp(u==units,statei,statej,:) = sq(mean(permpool(perminds)));
                else 
                    stsCor(u==units,statei,statej,:,:,2) = ones(2,2);
                end
                
            end
        end
    end

    dstruct(ses).clu   = pfs{1}.data.clu(:);
    dstruct(ses).elClu = pfs{1}.data.elClu(:);
    dstruct(ses).el    = pfs{1}.data.el(:);
    dstruct(ses).ratecorr = stsCor;   
    dstruct(ses).permsamp = permsamp; 
    dstruct(ses).overlap  = overlap;  
end
