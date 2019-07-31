
% ACCUMULATE phase precession stats ----------------------------------------------------------------
overwrite = false;
numIter = 100;
parmDRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsDRZall = zeros([unitCnt,2,numIter+1,numStates]);
rDRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoDRZall      = zeros([unitCnt,1,numIter+1,numStates]);
parmHRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsHRZall = zeros([unitCnt,2,numIter+1,numStates]);
rHRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoHRZall      = zeros([unitCnt,1,numIter+1,numStates]);

for t = 1:numel(Trials),

    Trial = Trials{t};    
    unitSubsetTempTemp = units{t};
    fprintf('Processing Trial: %s\nComputational time: ',Trial.filebase)        
   

% LOAD marker position from motion capture data
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);
    
% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
    Trial.lfp.filename = [Trial.name,'.lfp'];
    try, lfp = Trial.load('lfp',sessionList(t).thetaRef);
    catch, lfp = Trial.load('lfp',sessionList(t).thetaRef);
    end
    phz = lfp.phase([5,13]);    
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);



% COMPUTE placefield centered rate scaled distance metric with approach/depart signature
    ghz = compute_ghz(Trial,unitSubsetTemp,pfts{t},'sampleRate',sampleRate,'sigma',sigma);
    ddz = compute_ddz(Trial,unitSubsetTemp,pfts{t},'sampleRate',sampleRate);
    drz = compute_drz(Trial,unitSubsetTemp,pfts{t},'sampleRate',sampleRate);

% $$$ % COMPUTE head frame of reference vectors for all time points
% $$$     fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
% $$$     hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
% $$$     hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
% $$$     hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
% $$$ 
% $$$     tvec = circshift(fxyz(:,'hcom',[1,2]),-10)-circshift(fxyz(:,'hcom',[1,2]),10);
% $$$     tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
% $$$     tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);
% $$$ 
% $$$ % NAN sample points where spikes occured lateral to the head
% $$$     pfhr = nan([size(xyz,1),numel(unitSubsetTemp),2]);
% $$$     for u = 1:numel(unitSubsetTemp),%&pfts{t}.data.spar>0.15&pfts{t}.data.spar<0.3),
% $$$         [mxr,mxp] = pfts{t}.maxRate(unitSubsetTemp(u));
% $$$         pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
% $$$         hrz(abs(pfhr(:,u,2))>100,u) = nan;
% $$$     end
    
% ISOLATE approaching and retreating trajectories
% $$$     drzp = drz;    drzp(drzp<0)=nan;
% $$$     ddzp = ddz;    ddzp(ddzp<0)=nan;    
% $$$     drzn = drz;    drzn(drzn>0)=nan;
% $$$     ddzn = ddz;    ddzn(ddzn>0)=nan;    
% $$$     hrzp = hrz;    hrzp(hrzp<0)=nan;
% $$$     hdzp = ddz;    hdzp(hdzp<0)=nan;    
% $$$     hrzn = hrz;    hrzn(hrzn>0)=nan;
% $$$     hdzn = ddz;    hdzn(hdzn>0)=nan;        

% LOCATE unit row in group stats matrix
    if t == 1,  ind = 1:sessionUnitCnt(1);
    else,       ind = (sum(sessionUnitCnt(1:(t-1)))+1):sum(sessionUnitCnt(1:t));
    end

% COMPUTE phase precession stats for all states    
    for s = 1:numStates,
        spkpp = Trial.spk.copy();
        spkpp.create(Trial,xyz.sampleRate,states{s},unitSubsetTemp,'deburst');

% $$$         [parmDRZall(ind,:,:,s), phzStatsDRZall(ind,:,:,s), rDRZall(ind,:,:,s), rhoDRZall(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,drz,ddz,phz,spkpp,unitSubsetTemp,[],[],numIter,...
% $$$                                       ['DRZ-',states{s},'-',num2str(numIter)],overwrite);
        [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
            MjgER2016_phasePrecession(Trial,ghz,ddz,phz,spkpp,unitSubsetTemp,[],[],numIter,...
                                      ['GHZall-',states{s},'-',num2str(numIter)],overwrite);
% $$$         
% $$$         [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrz,ddz,phz,spkpp,unitSubsetTemp,[],[],numIter,...
% $$$                                       ['HRZall-',states{s},'-',num2str(numIter)],overwrite);
% $$$         [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrz,ddz,phz,spkpp,unitSubsetTemp,[],[],numIter,...
% $$$                                       ['HRZ-',states{s},'-',num2str(numIter)],overwrite);        

% HELD for later analysis
% $$$         [parmDRZpos(ind,:,s), phzStatsDRZpos(ind,:,s), rDRZpos(ind,:,s), rhoDRZpos(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spkpp,unitSubsetTemp);
% $$$         [parmDRZneg(ind,:,s), phzStatsDRZneg(ind,:,s), rDRZneg(ind,:,s), rhoDRZneg(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spkpp,unitSubsetTemp);
% $$$         [parmHRZpos(ind,:,s), phzStatsHRZpos(ind,:,s), rHRZpos(ind,:,s), rhoHRZpos(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzp,hdzp,phz,spkpp,unitSubsetTemp);
% $$$         [parmHRZneg(ind,:,s), phzStatsHRZneg(ind,:,s), rHRZneg(ind,:,s), rhoHRZneg(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzn,hdzn,phz,spkpp,unitSubsetTemp);
    end
        
    
end


% HELD for later analysis
% $$$ overwrite = false;
% $$$ numIter = 100;
% $$$ 
% $$$ parmDRZall     = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ phzStatsDRZall = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ rDRZall        = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ rhoDRZall      = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ parmHRZall     = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ phzStatsHRZall = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ rHRZall        = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ rhoHRZall      = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ parmHRZtrm     = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ phzStatsHRZtrm = zeros([unitCnt,2,numIter+1,numStates]);
% $$$ rHRZtrm        = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ rhoHRZtrm      = zeros([unitCnt,1,numIter+1,numStates]);
% $$$ 
% $$$ for t = 1:numel(Trials),
% $$$     Trial = Trials{t};    
% $$$     unitSubsetTemp = units{t};
% $$$     fprintf('Processing Trial: %s\nComputational time: ',Trial.filebase)        
% $$$     if t == 1,
% $$$         ind = 1:sessionUnitCnt(1);
% $$$     else,
% $$$         ind = (sum(sessionUnitCnt(1:(t-1)))+1):sum(sessionUnitCnt(1:t));
% $$$     end
% $$$ 
% $$$     for s = 1:numStates,
% $$$         [parmDRZall(ind,:,:,s), phzStatsDRZall(ind,:,:,s), rDRZall(ind,:,:,s), rhoDRZall(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,[],[],[],[],unitSubsetTemp,[],[],numIter,...
% $$$                                       ['DRZ-',states{s},'-',num2str(numIter)],overwrite);
% $$$         [parmHRZtrm(ind,:,:,s), phzStatsHRZtrm(ind,:,:,s), rHRZtrm(ind,:,:,s), rhoHRZtrm(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,[],[],[],[],unitSubsetTemp,[],[],numIter,...
% $$$                                       ['HRZ-',states{s},'-',num2str(numIter)],overwrite);
% $$$         [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,[],[],[],[],unitSubsetTemp,[],[],numIter,...
% $$$                                       ['HRZall-',states{s},'-',num2str(numIter)],overwrite);        
% $$$     end
% $$$ end


% END ACCUMULATE phase precession stats ------------------------------------------------------------
