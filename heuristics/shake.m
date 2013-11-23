function [shake_state,shake_feature] = shake(Session,method,varargin)
switch method
  case 'freq'
    [shakeThresh,nFFT,windowsize,display] = DefaultArgs(varargin,{0.0001,2^8,2^6,0});

    %% Temp Default Args
    %shakeThresh =0.0001;
    %nFFT =2^8;
    %windowsize = 2^6;

    %% Non-spectral Feature
    slpv = Session.ang(:,Session.Model.gmi('spine_lower'),Session.Model.gmi('pelvis_root'),1);
    slsm = Session.ang(:,Session.Model.gmi('spine_lower'),Session.Model.gmi('spine_middle'),1);

    slpvtp = slpv(1:end-1)-slpv(2:end);
    turns = zeros(size(slpvtp));
    turns(slpvtp>pi) = 1;
    turns(slpvtp<-pi) = -1;
    tcount = cumsum(turns);
    stcp = [0;tcount];
    cslpv = slpv+stcp*2*pi;

    slsmtp = slsm(1:end-1)-slsm(2:end);
    turns = zeros(size(slsmtp));
    turns(slsmtp>pi) = 1;
    turns(slsmtp<-pi) = -1;
    tcount = cumsum(turns);
    stcp = [0;tcount];
    cslsm = slsm+stcp*2*pi;

    shake = cslsm-cslpv;
    shake(isnan(shake))=0;
    wshake = WhitenSignal(shake);
    [y,f,t] = mtchglong(wshake,nFFT,Session.xyzSampleRate,windowsize,round(0.75*windowsize),[],[],[],[5 50]);
    shakep = mean(y(:,14<=f&f<=16),2);




    shake_feature = zeros(size(Session.xyz,1),1);
    tfet = d2t(shake_feature,Session.xyzSampleRate,0);
    [~,indconv] = min(abs(tfet-t(2)));
    shake_pre_feature = reshape(repmat(shakep',indconv,1),1,[]);
    while length(shake_pre_feature)>size(Session.xyz,1),
        indconv = indconv-1;
        shake_pre_feature = reshape(repmat(shakep',indconv,1),1,[]);
    end
    ignore_ind = [];
    while length(ignore_ind)+length(shake_pre_feature)~=length(shake_feature),
        fet_shift = length(shake_feature(shake_feature~=-1))-length(shake_pre_feature);
        if isprime(fet_shift), 
            ignore_ind = cat(1,ignore_ind,length(shake_feature)-fet_shift:length(shake_feature)-1);
        end
        ignore_ind = cat(1,ignore_ind,cumsum(round(size(Session.xyz,1)/fet_shift).*ones(fet_shift,1)));
        ignore_ind = unique(ignore_ind);
        ignore_ind = ignore_ind(ignore_ind<=length(shake_feature)); 
        shake_feature(ignore_ind) = -1;
    end
    shake_feature(shake_feature~=-1) = shake_pre_feature;
    if shake_feature(end)==-1, 
        shake_feature(end) = shake_feature(find(shake_feature~=-1,1,'last'));
    end    
    shake_feature(shake_feature==-1) = shake_feature(circshift(shake_feature,-1)==-1);


    shake_state = zeros(size(Session.xyz,1),1);
    shake_state(shake_feature>shakeThresh) = 1;
    shake_state = ThreshCross(shake_state,0.5,0);
end
