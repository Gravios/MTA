% [Fets PCs] = Feature(Spk, nChannels, SpkSampls, Spikes2Load, nPCA, OnePCASet)
%
% a matlab version of sfeature.
%
% Spk may be an array of spikes or a file name.
%
% nChannels will default to 4
% SpkSampls will default to 32
% Spikes2Load will default to inf - i.e. load them all.
% nPCA will default to 3.
% OnePCASet will do it so it uses the same principal components for all electrodes
%    -> default for this is 0, so it's like the usual sfeature.
%
% returns Fets - which is just like the feature file, and
% PCs, which is the principal component waveforms
% PCs will be SpkSampls by nPCA if OnePCASet is 1
% or SpkSampls by nPCA by nChannels if OnePCASet is 0
%
% Author: Ken Harris ( I think )
%

function [Fets, PCs] = fet_spk(Spk, nChannels, SpkSampls, Spikes2Load, nPCA, OnePCASet)


if (~isnumeric(Spk))
    if (nargin<2) nChannels = 4; end;
    if (nargin<3) SpkSampls = 32; end;
    if (nargin<4) Spikes2Load = inf; end;
	Spk = LoadSpk(Spk, nChannels, SpkSampls, Spikes2Load);
else
    nChannels = size(Spk, 1);
    SpkSampls = size(Spk, 2);
end;

if (nargin<5) nPCA = 3; end
if (nargin<6) OnePCASet = 0; end;

nSpikes = size(Spk,3);

if OnePCASet
	% make a single vector of all spikes
	fprintf('Rearranging data...\n');
	AllSpikes = reshape(permute(Spk,[3 1 2]), nSpikes*nChannels, SpkSampls);
	
	%subtract out mean
	AllSpikes = AllSpikes - repmat(mean(AllSpikes,1), nSpikes*nChannels, 1);
	
	% do the pca (i.e. svd)
	fprintf('Making covariance matrix...\n');
	CovMat = cov(AllSpikes);
	
	fprintf('Calculating eigenvectors...\n');
	[PCs, d] = eigs(CovMat, eye(size(CovMat)),nPCA);
	
	%now make feature matrix
	fprintf('Calculating features...\n');
	Fets = AllSpikes * PCs;
	Fets = reshape(Fets, [nSpikes,nChannels, nPCA]);
	Fets = permute(Fets, [1 3 2]);
	Fets = reshape(Fets, [nSpikes, nPCA*nChannels]);
else
	% go through channels
	for Channel=1:nChannels
		% calculate pca for this channel
		ChannelSpikes = squeeze(Spk(Channel, :, :))';
		CovMat = cov(ChannelSpikes);
        
        % calculate eigenvalues: why does it print by default?
        opts = struct('disp', 0);
		[ChannelPCs, d] = eigs(CovMat, eye(size(CovMat)),nPCA, 'LA', opts);
		
		% make it so PCs all have the same sign (hopefully...)
		ChannelPCs = ChannelPCs .* repmat(sign(ChannelPCs(13,:)), SpkSampls, 1);
		
		% now make feature vector
		ChannelFets = ChannelSpikes * ChannelPCs;
		
		% now add both of these to the output variables
		PCs(:,:,Channel) = ChannelPCs;
		Fets(:,(Channel-1)*nPCA + 1 : Channel*nPCA) = ChannelFets;
		
	end;
end;