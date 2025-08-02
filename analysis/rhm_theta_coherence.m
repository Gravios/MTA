configure_default_args();
MjgER2016_load_data();


tind = 20;
Trial = Trials{tind};
sampleRate = 250;
stc = copy(Trial.stc);
xyz = preproc_xyz(Trial,'trb',sampleRate);
vxy = vel(filter(copy(xyz),'ButFilter',4,2.4,'low'),{'spine_lower','nose','hcom'},[1,2]);
rhm = fet_rhm(Trial,sampleRate);
rhmPhz = phase(rhm,[5,12]);


vind = {};
vind{1} = true(size(mvxy(:,2)));
vind{2} = mvxy(:,2)<10&mvxy(:,2)>2;
vind{3} = mvxy(:,2)>10;


tlfp = Trial.load('lfp',sessionList(tind).subject.channelGroup.theta);
tlfp.resample(rhm);
cshifts = [0,-6500:125:6500];

for iter = 1:numel(cshifts)
    disp(iter)
    tic
    rhmCs = copy(rhm);
    rhmCs.data(nniz(rhmCs.data)) = circshift(rhm.data(nniz(rhmCs.data)),cshifts(iter));
    elfp = copy(tlfp);
    elfp.data = cat(2,elfp.data,rhmCs.data);
    specArgsTheta = struct('nFFT',2^9,...
                           'Fs',  elfp.sampleRate,...
                           'WinLength',2^8,...
                           'nOverlap',2^8*0.5,...
                           'NW',3,...
                           'Detrend',[],...
                           'nTapers',[],...
                           'FreqRange',[1,20]);
    [mys,mfs,mts] = fet_spec(Trial,elfp,'mtcsdglong',false,[],specArgsTheta,[],true);
    
    for sts = 1:numel(states)
        sind = logical(get(resample(cast([stc{states{sts}}],'TimeSeries'),mys),'data'));
        for vts = 1:numel(vind)
            ind = vind{vts}&sind;
            scoher(:,sts,vts,iter) = sum(abs(mys(ind,:,1,2)))./sum(sqrt(mys(ind,:,1,1).*mys(ind,:,2,2)));
        end
    end
    toc
end

