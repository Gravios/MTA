function [P,phzStats,Rmax,rho] = MjgER2016_phasePrecession(Trial,varargin)
%function [P,phzStats,Rmax] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units)
%
% Compute phase precession coefficients between distance restricted drz and lfp phase
% 
% Varargin:
%
%    drz:        directional rate zones
%    ddz:        directional distance zones
%    phz:        Local field potential phase
%    spk:        MTASpk object which holds spike events
%    units:      list of units for computation
%    distThresh: limit phase precession analysis to specified radius around placefield center
%    mResults:   number of best fits to be returned
%
% Output:
%
%    P - NumericMatrix[(number of units), 1, : 
%    phzStats - 
%    Rmax -
%    drzHCnt - 


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('drz',                            [],                                           ...
                 'ddz',                            [],                                           ...
                 'phz',                            [],                                           ...
                 'spk',                            [],                                           ...
                 'units',                          [],                                           ...
                 'distThresh',                     300,                                          ...
                 'fitRange',                       -pi:0.001:pi,                                 ...
                 'numIter',                        100,                                          ...
                 'tag',                            '',                                           ...
                 'overwrite',                      false                                         ...
);
[drz,ddz,phz,spk,units,distThresh,fitRange,numIter,tag,overwrite] =                              ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
if isempty(tag),
    tag = DataHash({spk,units,distThresh,fitRange,numIter});
end
%---------------------------------------------------------------------------------------------------

dataFilePath = fullfile(Trial.spath,[Trial.filebase,'-phasePrecession-',tag,'.mat']);



if exist(dataFilePath,'file') && ~overwrite,
    load(dataFilePath);
else,
    if isempty(phz),
        error('MTA_ANALYSIS_ERROR','MTA:analysis:MjgER2016:MjgER2016_phasePrecession:PhaseAbsent');
        %phz = phase(resample(Trial.load('lfp',[Trial.ephys.electrodes(:).CA1pyrThetaRef]),drz),[6,12]);
    end
    if isempty(drz),
        drz = compute_drz(Trial,units);
    end
    

    saveVars = {'P','phzStats','Rmax','rho'};

    % use unit 33 from Ed10-20140817.cof.gnd

    P          = nan([numel(units),2,numIter+1]);
    phzStats   = nan([numel(units),2,numIter+1]);
    Rmax       = nan([numel(units),1,numIter+1]);
    rho        = nan([numel(units),1,numIter+1]);

    %drzHCnt    = zeros([numel(units),10]);

    % GET res

    for uind = 1:numel(units);
        % GET unit index
        tic
        unit = units(uind);
        
        if ~isempty(spk),
            % GET spikes times of unit    
            % REMOVE spikes outside of position acquisition periods
            % REMOVE spikes outside of distance threshold
            res = spk(unit);
            res(res>size(drz,1))=[];
            res(drz(res,uind)==0)=[];            
            res (abs(ddz(res, uind))>=distThresh) = [];            
            % SKIP fitting parameters if too few spikes
            if numel(res)<10,
                continue;
            end        
            % GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
            % GET phase values at times of spikes   
            % IGNORE spikes where drz or phase are nans
            drzspk = drz(res,uind);
            phzspk = phz(res,1);            
        else,% compute phase precession for single unit
            drzspk = drz(abs(ddz)<=distThresh);
            phzspk = phz(abs(ddz)<=distThresh);
        end
        
        gind = ~isnan(drzspk)&~isnan(phzspk);

        % SKIP fitting parameters if too few spikes
        if sum(gind)<10,
            continue;
        end    
        
        %drzHCnt(uind,:) = histcounts(drzspk,linspace(-1,1,11));


        lin  = cat(3,drzspk(gind),permute(Shuffle(repmat(drzspk(gind),[1,numIter]),1),[1,3,2])); 
        circ = repmat(phzspk(gind),[1,1,numIter+1]);

        cosPart = [];
        try,  cosPart = sum(cos(bsxfun(@minus,circ,2*pi*bsxfun(@times,repmat(fitRange,[1,1,numIter+1]),lin))),1);
        catch err, disp(err);
            for iter = 1:numIter+1,
                cosPart(:,:,iter) = sum(cos(bsxfun(@minus,circ(:,1,iter),2*pi*bsxfun(@times,fitRange,lin(:,1,iter)))),1);
            end
        end

        sinPart = [];
        try,  sinPart = sum(sin(bsxfun(@minus,circ,2*pi*bsxfun(@times,repmat(fitRange,[1,1,numIter+1]),lin))),1);
        catch err, disp(err);
            for iter = 1:numIter+1,
                sinPart(:,:,iter) = sum(sin(bsxfun(@minus,circ(:,1,iter),2*pi*bsxfun(@times,fitRange,lin(:,1,iter)))),1);
            end
        end


        R = sqrt((cosPart./length(circ)).^2+...
                 (sinPart./length(circ)).^2 );

        %[lmi,~] = LocalMinima(-R',0,0,mResults);
        %lmid = find(~isnan(lmi));
        
        [~,lmi] = min(-R,[],2);

        % Rmax fit quality 
        %Rmax(uind,1,lmid) = R(lmi(lmid));
        lmind = sq(lmi)+(0:size(R,3)-1)'.*size(R,2);
        Rmax(uind,1,:) = R(lmind);  
        % P:  Regression Parm: [Slope,Offset]
        P(uind,:,:) = permute([2*pi*fitRange(lmi);atan2(sinPart(lmind),cosPart(lmind))'],[3,1,2]);

        linC = mod(multiprod(abs(P(uind,1,:)),lin,[1],[1]),2*pi);
        circMean = atan2(sum(sin(circ)),sum(cos(circ)));
        linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
        rho(uind,1,:) = sum(sin(circ-circMean).*sin(linC-linCMean))...
            ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
        

% $$$ clf
% $$$ subplot(231);
% $$$ hist(sq(P(uind,1,2:end)),-10:0.1:10),Lines(P(uind,1,1),[],'r');   axis('tight');
% $$$ subplot(232);
% $$$ hist(sq(rho(uind,1,2:end)),-1:0.01:1),Lines(rho(uind,1,1),[],'r');   axis('tight');
% $$$ subplot(233);
% $$$ hist(sq(Rmax(uind,1,2:end)),0:0.001:1),Lines(Rmax(uind,1,1),[],'r');axis('tight');
% $$$ subplot(235);hold('on');
% $$$ plot(lin(:,1,2),phzspk(gind),'.b')
% $$$ plot(lin(:,1,2),phzspk(gind)+2*pi,'.b')
% $$$ subplot(236);hold('on');
% $$$ plot(drzspk(gind),phzspk(gind),'.b')
% $$$ plot(drzspk(gind),phzspk(gind)+2*pi,'.b')
% $$$ waitforbuttonpress();
        % Collect residuals of the theta model for each state
% $$$     phi = nan([numel(lin),mResults]);
% $$$     for r = lmid',
% $$$         [phi(:,r)] = polyval([2*pi;1].*sq(P(uind,:,:))',lin);
% $$$     end
% $$$     residuals = bsxfun(@minus,circ,phi);
% $$$     residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
% $$$     residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;
% $$$ 
% $$$     phzStats(uind,:,:) = [circ_mean(residuals);circ_std(residuals)];
        
        % COMPUTE circular-linear correlation coefficient
% $$$     linC = mod(abs(P(uind,1,1))*lin,2*pi);
% $$$     circMean = atan2(sum(sin(circ)),sum(cos(circ)));
% $$$     linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
% $$$     rho(uind,1,lmid) = sum(sin(circ-circMean).*sin(linC-linCMean))...
% $$$         ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));

        toc
    end

    save(dataFilePath,saveVars{:},'-v7.3');
end
% $$$ figure,
% $$$ hold('on');
% $$$ plot([lin,lin],[circ,circ+2*pi],'b.')
% $$$ plot([-1,1],P(uind,1,1)*[-1,1]+P(uind,2,1),'-m','LineWidth',1)
% $$$ plot([-1,1],P(uind,1,1)*[-1,1]+P(uind,2,1)+2*pi,'-m','LineWidth',1)
