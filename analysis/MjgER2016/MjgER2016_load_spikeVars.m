function [spkv,statesSpk] = MjgER2016_load_spikeVars(Trials,units,sessionList,varargin)
% MjgER2016_load_spikeVars
%  Status: perm
%  Type: Analysis
%  Author: Justin Graboski
%  Description: group phase precession stats
%

%states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('tag',                 '',                                                      ...
                 'statesSpk',           {{'theta','loc','pause','rear',                          ...
                                          'hloc','hpause','lloc','lpause',                       ...
                                          'groom','sit'}},                                       ...
                 'sampleRate',          250,                                                     ...
                 'sigma',               150,                                                     ...
                 'overwrite',           false                                                    ...
);
[tag,statesSpk,sampleRate,sigma,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash({statesSpk,sigma,sampleRate});
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
filename = fullfile(MTASession([]).path.project,'analysis',[mfilename,'-',tag,'.mat']);


if ~exist(filename,'file') || overwrite,
    
    spkmap = [];
    spkdrz = [];
    spkddz = [];    
    spkhrz = [];
    spkghz = [];
    spkpch = [];
    spkego = [];
    spkphz = [];
    spkvxy = [];
    spkstc = [];
    spkfbr = [];
    spksvd = [];
    spkavl = [];
    spktrans = repmat({[]},[2,numel(statesSpk)]);

    for tind = 1:numel(Trials)
        Trial = Trials{tind};
        unitSet = units{tind};

% LOAD xyz trajectories
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('ButFilter',3,30,'low');    
        xyz.resample(sampleRate);        

% COMPUTE speed of body and head
        vxy = vel(filter(copy(xyz),'ButFilter',3,30,'low'),{'spine_lower','nose','hcom'},[1,2]);

% COMPUTE body and head pitch
        fet = fet_HB_pitchB(Trial,sampleRate);

% COMPUTE angular velocity
        bxyz = xyz.copy();
        bxyz.data = bxyz(:,{'spine_lower','spine_upper','head_back','head_front'},:);
        bxyz.model = bxyz.model.rb({'spine_lower','spine_upper','head_back','head_front'});
        ang = create(MTADang,Trial,bxyz);
        avl = MTADfet.encapsulate(Trial,...
                                  [circ_dist(circshift(ang(:,3,4,1),-3),...
                                             circshift(ang(:,3,4,1),3)).*xyz.sampleRate,...
                            circ_dist(circshift(ang(:,1,2,1),-3),...
                                      circshift(ang(:,1,2,1),3)).*xyz.sampleRate],...
                                  ang.sampleRate,...
                                  'angular velocity',...
                                  'angvel',...
                                  'v');
        
% LOAD behavioral states
        stc = Trial.stc.copy();
        
% LOAD units
        spk = Trial.load('spk',sampleRate,'gper',unitSet,'deburst');
        pft = pfs_2d_theta(Trial,unitSet);
        hrz = compute_hrz(Trial,unitSet,pft,'sampleRate',sampleRate);
        ghz = compute_ghz(Trial,unitSet,pft,'sampleRate',sampleRate,'sigma',sigma);
        ddz = compute_ddz(Trial,unitSet,pft,'sampleRate',sampleRate);
        drz = compute_drz(Trial,unitSet,pft,'sampleRate',sampleRate);
        
% LOAD phase
        Trial.lfp.filename = [Trial.name,'.lfp'];
        try,  lfp = Trial.load('lfp',sessionList(tind).thetaRefGeneral);
        catch,lfp = Trial.load('lfp',sessionList(tind).thetaRefGeneral);
        end    
        phz = lfp.phase([5,13]);    
        phz.data = unwrap(phz.data);
        phz.resample(xyz);    
        phz.data = mod(phz.data+pi,2*pi)-pi;
        lfp.resample(xyz);

% COMPUTE head vectior        
        hvec = xyz(:,'head_front',[1,2])-xyz(:,'head_back',[1,2]);
        hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
        hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);        
        
% CONVERT state periods to state matrix
        stcm = stc2mat(stc,xyz,statesSpk);
        
% COLLATE state transitions
        ststrans = {};
        for state = 1:numel(statesSpk),
            samples = round(xyz.sampleRate/2);
            stsper = [Trial.stc{statesSpk{state}}];
            stsdiff = [[stsper.data(:,2);size(xyz,1)]-[0;stsper.data(:,1)]]';
            for s = 1:2,
                trans = stsper.data(:,s)';
                trans(stsdiff(s:end-2+s)<=samples) = [];
                trans((trans-samples)<=0|(trans+samples)>=size(xyz,1)) = [];
                ststrans{s,state} = xyz.copy('empty');
                ststrans{s,state}.data = nan([size(xyz,1),1]);
                
                for transition = trans,
                    ststrans{s,state}.data(transition-samples:transition+samples) = ...
                        linspace(-1,1,2*samples+1);
                end
            end
        end

% COLLECT vars for all spikes        
        for unit = 1:numel(unitSet);
% GET unit index    
% GET spikes times of unit    
% REMOVE spikes outside of position acquisition periods
% REMOVE spikes outside of distance threshold
            [mxr,mxp] = pft.maxRate(unitSet(unit));
            
            res = spk(unitSet(unit));
            res(res>size(ddz,1))=[];
% GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
% GET phase values at times of spikes   
% IGNORE spikes where drz or phase are nans
            spkmap = cat(1,spkmap,[ones([numel(res),1]).*tind,...
                                ones([numel(res),1]).*unitSet(unit)]);
            
            spkego = cat(1,spkego,multiprod(bsxfun(@minus,mxp,sq(xyz(res,'hcom',[1,2]))),hvec(res,:,:),2,[2,3]));
            spkstc = cat(1,spkstc,stcm(res,:));
            spkpch = cat(1,spkpch, fet(res,:));        
            spkvxy = cat(1,spkvxy, vxy(res,:));        
            spkavl = cat(1,spkavl, avl(res,:));        
            spkdrz = cat(1,spkdrz, drz(res,unit));
            spkhrz = cat(1,spkhrz, hrz(res,unit));
            spkddz = cat(1,spkddz, ddz(res,unit));
            spkghz = cat(1,spkghz, ghz(res,unit));
            spkphz = cat(1,spkphz, phz(res,1));        
            for state = 1:numel(statesSpk),
                for s = 1:2,
                    spktrans{s,state} = cat(1,spktrans{s,state},ststrans{s,state}(res));
                end
            end
        end% for unit
    end%for tind

    spkv.states     = statesSpk;
    spkv.sampleRate = sampleRate;
    spkv.sigma      = sigma;

    spkv.map   = spkmap;
    spkv.drz   = spkdrz;
    spkv.ddz   = spkddz;    
    spkv.hrz   = spkhrz;
    spkv.ghz   = spkghz;
    spkv.pch   = spkpch;
    spkv.ego   = spkego;
    spkv.phz   = spkphz;
    spkv.vxy   = spkvxy;
    spkv.stc   = spkstc;
    spkv.svd   = spksvd;
    spkv.avl   = spkavl;
    spkv.trans = spktrans;

    save(filename,'spkv','-v7.3');
else
    load(filename);
end
% END MAIN -----------------------------------------------------------------------------------------