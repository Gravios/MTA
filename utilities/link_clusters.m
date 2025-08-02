% LinkClusters(FileBase1, FileBase2)
%
% compares two sets of clustered data to see if the clusters
% correspond to the same cells. 
%
% It does this by looking at the waveforms (not the .fet files
% as these may be wrt a different pc basis), doing a pca (all
% electrodes at once) and then comparing with MANOVA.
%
% Author: Ken Harris ( I think )
%

function link_clusters(FileBase1, FileBase2)

Session = MTASession([]);
filebase1 = fullfile(Session.path.project,FileBase1,FileBase1);
filebase2 = fullfile(Session.path.project,FileBase2,FileBase2);

Par1 = LoadPar([filebase1 '.xml']);
Par2 = LoadPar([filebase2 '.xml']);

if Par1.nElecGps~=Par2.nElecGps
    warning('different numbers of tetrodes on the 2 files')
end

allSpkStats = {};

cluCnt = [];
for i=1:min(Par1.nElecGps, Par2.nElecGps)

    fprintf('\nElec %d', i);

    Clu1 = LoadClu([filebase1 '.clu.' num2str(i)]);
    Clu2 = LoadClu([filebase2 '.clu.' num2str(i)]);
    nClu1 = max(Clu1);
    nClu2 = max(Clu2);

    cluCnt = cat(2,cluCnt,[max(Clu1);max(Clu2)]);

    if min(nClu1, nClu2)==1
        fprintf(' No clusters!');
        continue
    end
    
% $$$     if ~isequal(rmfield(Par11, 'FileName'), rmfield(Par12, 'FileName'))
% $$$         warning('different .par. files');
% $$$     end
    
% $$$     fprintf(' Loading Spikes ');
% $$$     
% $$$     Spk1 = LoadSpk([filebase1 '.spk.' num2str(i)],   ... path to spk file
% $$$                     numel(Par1.SpkGrps(i).Channels), ... number of channels in spk group
% $$$                     Par1.SpkGrps(i).nSamples);         % number of samples per waveform
% $$$ 
% $$$     Spk2 = LoadSpk([filebase2 '.spk.' num2str(i)],   ... path to spk file
% $$$                     numel(Par2.SpkGrps(i).Channels), ... number of channels in spk group
% $$$                     Par2.SpkGrps(i).nSamples);         % number of samples per waveform
% $$$     
% $$$     Size1 = size(Spk1);
% $$$     Size2 = size(Spk2);
 
% $$$     AllSpk = double(cat(3, Spk1, Spk2));
% $$$     Fet = fet_spk(AllSpk);

    Fet1 = LoadFet([filebase1 '.fet.' num2str(i)]); % path to fet file
    Fet2 = LoadFet([filebase2 '.fet.' num2str(i)]); % path to fet file

    Fet = [Fet1(:,1:end-1);Fet2(:,1:end-1)];
    % Assign group Ids    
    gp = [zeros([numel(Clu1),1]) ; ones([numel(Clu2),1])];    

    
    fprintf('Processing ');
    clear d p stats
    for c1=2:nClu1
        for c2=2:nClu2
            % find inds for each clu and concatenate 
            ind = [Clu1==c1;Clu2==c2];            
            % Generate stats from spike features
            [d(c1,c2), p(c1,c2) stats(c1, c2)] = manova1(Fet(ind,:), gp(ind));
            % prove the computer is not being lazy 
            fprintf('.');
        end
    end
     
    allSpkStats{i} = reshape([stats.lambda], [nClu1-1 nClu2-1]);
    figure,imagesc(2:nClu1, 2:nClu2, allSpkStats{i}', [0 1]);
    colorbar;
    pause

    

    s = [4,15];
    figure,
    subplot(121),
    imagesc(reshape(dp(1).Pfs.data.rateMap(:,cluCnt(1)+s(1),1),cellfun(@numel,dp(1).Pfs.adata.bins))') 
    subplot(122),
    imagesc(reshape(dp(2).Pfs.data.rateMap(:,cluCnt(2)+s(2),1),cellfun(@numel,dp(2).Pfs.adata.bins))') 



end