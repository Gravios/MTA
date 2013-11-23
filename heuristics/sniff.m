function [sniff_state,sniff_feature] = sniff(Session,method,varargin)
[sniffThresh,nFFT,windowsize,freqRange,display] = DefaultArgs(varargin,{10^-3.5,2^10,2^8,[0.1,20],0});
switch method
  case 'freq'
    dsniff = sqrt(sum((Session.xyz(:,Session.Model.gmi('head_back'),:)-Session.xyz(:,Session.Model.gmi('spine_upper'),:)).^2,3));
    wdsniff = WhitenSignal(dsniff);

    [y,f,t] = mtchglong(wdsniff,1024,Session.xyzSampleRate,windowsize,round(0.75*windowsize),[],[],[],[0.1 20]);
    sniffp = mean(y(:,8<f&f<14),2);

% $$$     shiftedSync = Session.syncPeriods-Session.syncPeriods(1,1);
% $$$     spectralBoundaries = [ shiftedSync(1:end-1,2),shiftedSync(2:end,1)];
% $$$     [~,timeBoundaries] =  NearestNeighbour(t+windowsize/(2*Session.xyzSampleRate), spectralBoundaries/Session.lfpSampleRate,'both');
% $$$     timeBoundaries = [timeBoundaries([1:2:end-1]);timeBoundaries([2:2:end])];
% $$$     timeBoundaries = reshape(timeBoundaries,2,[])';
% $$$ 
% $$$     for i = 1:size(timeBoundaries,1)
% $$$         sniffp(timeBoundaries(i,1)-5:timeBoundaries(i,2)+5)=mean(sniffp);
% $$$     end

    sniff_feature = zeros(size(Session.xyz,1),1);

    tfet = d2t(sniff_feature,Session.xyzSampleRate,0);
    [~,indconv] = min(abs(tfet-t(2)));
    sniff_pre_feature = reshape(repmat(sniffp',indconv,1),1,[]);
    while length(sniff_pre_feature)>size(Session.xyz,1),  
        indconv = indconv-1;      
        sniff_pre_feature = reshape(repmat(sniffp',indconv,1),1,[]);
    end
    ignore_ind = [];
    while length(ignore_ind)+length(sniff_pre_feature)~=length(sniff_feature),
        fet_shift = length(sniff_feature(sniff_feature~=-1))-length(sniff_pre_feature);


        if isprime(fet_shift), 
            ignore_ind = cat(1,ignore_ind,(length(sniff_feature)-fet_shift:length(sniff_feature)-1)');
        end

        ignore_ind = cat(1,ignore_ind,cumsum(round(size(Session.xyz,1)/fet_shift).*ones(fet_shift,1)));
        ignore_ind = unique(ignore_ind);
        ignore_ind = ignore_ind(ignore_ind<=length(sniff_feature)); 
        sniff_feature(ignore_ind) = -1;
    end

    sniff_feature(sniff_feature~=-1) = sniff_pre_feature;
    if sniff_feature(end)==-1, 
        sniff_feature(end) = sniff_feature(find(sniff_feature~=-1,1,'last'));
    end    

    sniff_feature(sniff_feature==-1) = sniff_feature(circshift(sniff_feature,-1)==-1);


    sniff_state = zeros(size(sniffp));
    sniff_state(sniffp>sniffThresh) = 1;
    sniff_state(sniffp<sniffThresh) = 0;
    sniff_state = ThreshCross(sniff_state,0.5,0);
    
  otherwise
    error(['Method: ' method ' does not exist'])
end


