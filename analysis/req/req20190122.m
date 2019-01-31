% req20190122.m
% compute time shifted behavior field information 
% OUTPUT : bhvInfoTS
%
% fields were computed in req20180611_tshift.m

global MTA_PROJECT_PATH

% LOAD data ---------------------------------------------------------------------------%
if ~exist('overwrite','var'),  overwrite = false;  end

if ~exist('sessionList','var'),  
    MjgER2016_load_data();

    % SET analysis parameters
    sampleRate = 250;   % Hz
end
    
dataFilePath = fullfile(MTA_PROJECT_PATH,                                       ...
                        'analysis',                                             ...
                        ['req20190122-data-',                                   ...
                         DataHash({[sessionList.sessionName],sampleRate}),'.mat']);


% COMPILE results ---------------------------------------------------------------------%

if ~exist(dataFilePath) || overwrite,
    Trial = Trials{1}; 
    unitSubset = units{1};        
    pfs = MTAApfs(Trial,'tag',['hbpptbpFS1v3_ts',num2str(-250)]);

    tshifts = -250:50:250;
    numPhzBins = pfs.adata.binSizes(2);
    eSpi = zeros([0,pfs.adata.binSizes(2),numel(tshifts)]);

    unitCount = 1;
    for t = 1:23,
        Trial = Trials{t}; 
        unitSubset = units{t};        
        
        pfs = {};
        for shift = tshifts,
            pfs{end+1} = MTAApfs(Trial,'tag',['hbpptbpFS1v3_ts',num2str(shift)]);
        end
        
        for unit = unitSubset
            for s = 1:numel(tshifts),
                rmap = plot(pfs{s},unit,1,'',[],false);
                rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=numPhzBins;
                zmap = reshape(permute(rmap,[1,3,2]),[],numPhzBins);
                zmap(repmat(~rmapNind(:),[1,numPhzBins])) = 0;            
                zmap(zmap(:)<0) = 0;
                rmapNind = rmapNind&reshape(sum(reshape(zmap,[],numPhzBins)','omitnan')~=0,...
                                            size(rmap,1),size(rmap,3));
                for p = 1:numPhzBins,
                    mrate = mean(reshape(zmap(:,p),[],1),'omitnan');
                    pmap = zmap(:,p);
                    pmap = pmap(rmapNind);
                    bhvInfoTS(unitCount,p,s) = sum(pmap./mean(pmap,'omitnan')...
                                              .*log2(pmap./mean(pmap,'omitnan')),'omitnan')./numel(pmap);
                end
            end   
            unitCount = unitCount+1;
        end
    end
    save(dataFilePath,'bhvInfoTS','tshifts');
else
    load(dataFilePath);
end




