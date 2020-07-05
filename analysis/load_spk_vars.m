function [spkv] = load_spk_vars(sessionList,Trials,units,varargin)
% MjgER2016_load_spikeVars
%  Status: perm
%  Type: Analysis
%  Author: Justin Graboski
%  Description: accumulate vars at spike times
%  
%  map - 

%states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('tag',                 '',                                                      ...
                 'states',           {{'theta','loc','pause','rear',                             ... 
                                         'hloc','hpause','lloc','lpause',                        ...
                                          'groom','sit'}},                                       ...
                 'rotation',                 zeros(size(sessionList)),                                ...
                 'sampleRate',          250,                                                     ...
                 'sigma',               150,                                                     ...
                 'overwrite',           false                                                    ...
);
[tag,states,rotation,sampleRate,sigma,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash({states,sigma,sampleRate});
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
filename = fullfile(MTASession([]).path.project,'analysis',[mfilename,'-',tag,'.mat']);

if ~exist(filename,'file') || overwrite,
    
    spkmap = []; % spike map relating clus to session id
    
    spkdrz = []; % Directed rate zones(data)
    spkddz = []; % distance directed zones (linear)
    spkhrz = []; % head-directed rate zones (data)
    spkghz = []; % gausian head-directed zones (gaussian)
    
    spkego = []; % position of place field center in head frame of reference
    spkbdy = [];
    spktrj = []; % position of place field center in trajectory frame of reference    

    spkhba = []; % yaw between head and body vectors
    spkhva = []; % angular velocity of body and head    
    spkhvl = [];
    spkhvf = [];    
    spkhbp = []; % pitch of head and body
    spkhbd = []; % distance between head and body
    spkvxy = []; % xy speed of the body and head
    spkbma = []; % body maze angle 
    
    spkphz = []; % theta phase
    spkstc = []; % behavioral state 

    for tind = 1:numel(Trials)
        Trial = Trials{tind};
        unitSet = units{tind};
        rot = rotation(tind);

% LOAD xyz trajectories
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('ButFilter',3,30,'low');    
        xyz.resample(sampleRate);        

% COMPUTE speed of body and head
        vxy = vel(xyz,{'hcom','spine_lower'},[1,2]);

% COMPUTE body and head pitch
        fet = fet_HB_pitchB(Trial,sampleRate);

% HBA - head-body angle
        hba = filter(copy(xyz),'ButFilter',3,30,'low');    
        xycoor = cat(2,...
                     hba(:,'spine_upper',[1,2])-hba(:,'bcom',[1,2]),...
                     hba(:,'nose',[1,2])-hba(:,'hcom',[1,2]));
        hba.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
        hba.data = circ_dist(hba.data(:,2),hba.data(:,1));
      
% HVA head angular velocity
        hva = filter(copy(xyz),'ButFilter',4,2,'low');
        xycoor = cat(2,...
                     hva(:,'spine_upper',[1,2])-hva(:,'bcom',[1,2]),...
                     hva(:,'nose',[1,2])-hva(:,'hcom',[1,2]));
        hva.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
        % Positive: CCW (Left)     Negative: CW (Right)
        hva.data = circ_dist(circshift(hva.data(:,2),-10),...
                               circshift(hva.data(:,2),+10));
        
% HBD - head-body distance
        hbd = copy(xyz);
        hbd.data = sqrt(sum(diff(hbd(:,{'hcom','bcom'},[1,2]),1,2).^2,3));

% HVL - head lateral movement
        %hvfl = fet_href_HXY(Trial,sampleRate,false,'trb',4);
        hvfl = fet_href_HXY(Trial,sampleRate,false,'trb',4,rot);
        
        
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
        hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
        hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
        hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
        hvec = multiprod(hvec,...
                         [cos(rot),-sin(rot);sin(rot),cos(rot)],...
                         [2,3],...
                         [1,2]);

        bvec = xyz(:,'spine_upper',[1,2])-xyz(:,'bcom',[1,2]);
        bvec = sq(bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3))));
        bvec = cat(3,bvec,sq(bvec)*[0,-1;1,0]);
% $$$         bvec = multiprod(bvec,...
% $$$                          [cos(rot),-sin(rot);sin(rot),cos(rot)],...
% $$$                          [2,3],...
% $$$                          [1,2]);
        
        tvec = circshift(xyz(:,'hcom',[1,2]),-round(sampleRate*0.1)) - ...
               circshift(xyz(:,'hcom',[1,2]), round(sampleRate*0.1));
        tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
        tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);
        
% COMPUTE 
        bodyMazeAng = [-xyz(:,'bcom',[1,2]),xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2])];
        bodyMazeAng = sq(bsxfun(@rdivide,bodyMazeAng,sqrt(sum(bodyMazeAng.^2,3))));
        bodyMazeAng = cart2pol(bodyMazeAng(:,:,1),bodyMazeAng(:,:,2));
        bodyMazeAng = circ_dist(bodyMazeAng(:,1),bodyMazeAng(:,2));

% CONVERT state periods to state matrix
        stcm = stc2mat(stc,xyz,states);
        
% $$$ % COLLATE state transitions
% $$$         ststrans = {};
% $$$         for state = 1:numel(states),
% $$$             samples = round(xyz.sampleRate/2);
% $$$             stsper = [Trial.stc{states{state}}];
% $$$             stsdiff = [[stsper.data(:,2);size(xyz,1)]-[0;stsper.data(:,1)]]';
% $$$             for s = 1:2,
% $$$                 trans = stsper.data(:,s)';
% $$$                 trans(stsdiff(s:end-2+s)<=samples) = [];
% $$$                 trans((trans-samples)<=0|(trans+samples)>=size(xyz,1)) = [];
% $$$                 ststrans{s,state} = xyz.copy('empty');
% $$$                 ststrans{s,state}.data = nan([size(xyz,1),1]);
% $$$                 
% $$$                 for transition = trans,
% $$$                     ststrans{s,state}.data(transition-samples:transition+samples) = ...
% $$$                         linspace(-1,1,2*samples+1);
% $$$                 end
% $$$             end
% $$$         end

% COLLECT vars for all spikes        
        res = [];
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
            
            spkego = cat(1,spkego,multiprod(bsxfun(@minus,...
                                                   mxp,...
                                                   sq(xyz(res,'hcom',[1,2]))),...
                                            hvec(res,:,:),2,[2,3]));
            spkbdy = cat(1,spkbdy,multiprod(bsxfun(@minus,...
                                                   mxp,...
                                                   sq(xyz(res,'hcom',[1,2]))),...
                                            bvec(res,:,:),2,[2,3]));
            spktrj = cat(1,spktrj,multiprod(bsxfun(@minus,...
                                                   mxp,...
                                                   sq(xyz(res,'hcom',[1,2]))),...
                                            tvec(res,:,:),2,[2,3]));
            spkstc = cat(1,spkstc,stcm(res,:));
            spkhbp = cat(1,spkhbp, fet(res,:));
            spkhva = cat(1,spkhva, hva(res,:));
            spkhba = cat(1,spkhba, hba(res,:));
            spkbma = cat(1,spkbma, bodyMazeAng(res,:));                
            spkhbd = cat(1,spkhbd, hbd(res,:));
            spkvxy = cat(1,spkvxy, vxy(res,:));
            spkhvf = cat(1,spkhvf, hvfl(res,1));
            spkhvl = cat(1,spkhvl, hvfl(res,2));
            spkdrz = cat(1,spkdrz, drz(res,unit));
            spkhrz = cat(1,spkhrz, hrz(res,unit));
            spkddz = cat(1,spkddz, ddz(res,unit));
            spkghz = cat(1,spkghz, ghz(res,unit));
            spkphz = cat(1,spkphz, phz(res,1));
        end% for unit
    end%for tind

    spkv.states     = states;
    spkv.sampleRate = sampleRate;
    spkv.sigma      = sigma;

    spkv.map   = spkmap;
    
    spkv.drz   = spkdrz;
    spkv.ddz   = spkddz;    
    spkv.hrz   = spkhrz;
    spkv.ghz   = spkghz;
    
    spkv.hva   = spkhva;
    spkv.hba   = spkhba;
    spkv.hvf   = spkhvf;    
    spkv.hvl   = spkhvl;    
    spkv.hbd   = spkhbd;    
    spkv.hbp   = spkhbp;
    spkv.vxy   = spkvxy;
    spkv.bma   = spkbma;
    
    spkv.ego   = spkego;
    spkv.trj   = spktrj;
    
    spkv.phz   = spkphz;
    spkv.stc   = spkstc;

    save(filename,'spkv','-v7.3');
else
    load(filename);
end
% END MAIN -----------------------------------------------------------------------------------------