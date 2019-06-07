%req20181220
% Compute idependent trajectory count within all regions of the behavioral ratemaps

global MTA_PROJECT_PATH;

MjgER2016_load_data();

% SET analysis parameters
sampleRate = 16;   % Hz

t = 20;    
Trial = Trials{t}; 
unitSubset = units{t};        

stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
hvec = xyz(:,'head_front',[1,2])-xyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);


pft = pfs_2d_theta(Trial,unitSubset);
[drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
[ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);

stcm = stc2mat(stc,xyz,states);

spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    

aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);

tbinInds = discretize(1:size(xyz,1),1:sampleRate*3:size(xyz,1));
for unit = unitSubset,
    trajCount = numel(nonzeros(diff(unique(tbinInds(aper & abs(drz(:,unit==unitSubset))<0.8 ...
                                      & abs(ddz(:,unit==unitSubset))<250)))-1));
end

