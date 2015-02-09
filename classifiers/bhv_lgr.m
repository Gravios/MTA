function Stc = bhv_lgr(Trial,Stc,varargin)
[train,states,model_filename,display] = DefaultArgs(varargin,{false,{},'MTA_manual_mknsrw_LGR_model.mat',true});


lrfet = fet_lgr(Trial);

if train   
    if isempty(states),   states = Trial.stc.list_state_attrib('label'); end

    Model_Information = struct(...
        'description',           '',                            ...
        'StcMode',               Trial.stc.mode,                ...
        'StcFilename',           Trial.stc.filename,            ...
        'SessionFilebase',       Trial.filebase,                ...
        'state_labels',          Trial.stc.list_state_attrib('label'),...
        'state_keys',            Trial.stc.list_state_attrib('key'));

    smat = max(stc2mat(Stc,ufet),[],2);
    ind = any(smat,2);
    B = mnrfit(lrfet(ind,:),smat(ind),'model','nominal');
    save(fullfile(fileparts(mfilename('fullpath')),model_filename),...
         'B','Model_Information');
else
    load(model_filename);
end    



d_state = mnrval(B,lrfet.data);

[~,w] = max(d_state,[],2);


[~,maxState] = max(d_state,[],2);

for i = 1:numel(states),
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,10),...
             Trial.xyz.sampleRate,...
             Trial.xyz.sync.copy,...
             Trial.xyz.origin,...
             states{i},...
             keys{i},...
             'TimePeriods');
end



% $$$ perind =[31600,33700];
% $$$ 
% $$$ figure,
% $$$ sp(1) = subplot(9,1,[1:4]);
% $$$ imagesc((1:size(ufet,1))/ufet.sampleRate,1:7,nunity(lrfet(:,1:7))');caxis([0,2]),colormap jet
% $$$ 
% $$$ % $$$ sp(2) = subplot(9,1,5);
% $$$ % $$$ rper = Trial.stc{'r'};
% $$$ % $$$ wper = Trial.stc{'w'};
% $$$ % $$$ rind(1) = find(rper(:,1)>perind(1)-1000,1,'first');
% $$$ % $$$ rind(2) = find(rper(:,2)<perind(2)+1000,1,'last');
% $$$ % $$$ wind(1) = find(wper(:,1)>perind(1)-1000,1,'first');
% $$$ % $$$ wind(2) = find(wper(:,2)<perind(2)+1000,1,'last');
% $$$ sp(2) = subplot(9,1,[6:9]);
% $$$ plot((1:size(ufet,1))/ufet.sampleRate,y)
% $$$ ylim([-.1,1.1])
% $$$ 
% $$$ hold on,Lines(Trial.stc{'w',1}(:),[],'k');
% $$$ hold on,Lines(Trial.stc{'r',1}(:),[],'r');
% $$$ hold on,Lines(Trial.stc{'m',1}(:),[],'m');
% $$$ 
% $$$ linkaxes(sp,'x');
% $$$ xlim(perind/xyz.sampleRate)
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ ind = logical(aper.data);
% $$$ % $$$ figure,
% $$$ % $$$ x = abs(circ_dist(dsa(ind,'spine_lower','spine_upper',1),dsa(ind,'spine_middle','head_front',1)));   xl = 0:.04:pi;
% $$$ % $$$ %x = dsa(ind,'spine_lower','head_front',3);   xl = 20:3:300;
% $$$ % $$$ y = dsx(ind,'spine_lower',3);   yl = 1:1:70;    
% $$$ % $$$ %y = dsa(ind,'spine_lower','pelvis_root',2);   yl = .2:.01:1.6;    
% $$$ % $$$ hist2([x,y],xl,yl),caxis([0,400])
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ figure,
% $$$ % $$$ for i = 1:3,
% $$$ % $$$     
% $$$ % $$$     x = circ_dist(dsa(w==i&aper,'spine_lower','spine_upper',1),dsa(w==i&aper,'spine_middle','head_front',1));   xl = 0:.05:pi;
% $$$ % $$$     %x = dsa(w==i&aper,'spine_lower','pelvis_root',2);   xl = .2:.01:1.6;
% $$$ % $$$     %x = dsa(w==i&aper,'pelvis_root','spine_middle',2);   xl = -.6:.01:1.4;
% $$$ % $$$     %x = rhm(w==i&aper.data);   xl = -1.5:.05:1.5;
% $$$ % $$$     %x = dsa(w==i&aper,'spine_lower','head_front',3);   xl = 20:3:300;
% $$$ % $$$     %x = ufet(w==i&aper);   xl = 18:.5:65;
% $$$ % $$$     y = dsa(w==i&aper,'spine_lower','pelvis_root',2);   yl = .2:.01:1.6;    
% $$$ % $$$     %y = dsx(w==i&aper,'spine_lower',3);                yl = 1:1:70;
% $$$ % $$$     %y = vel(w==i&aper,'head_front');                    yl = -1.5:.05:2;
% $$$ % $$$     subplotfit(i,4);
% $$$ % $$$     hist2([x,y],xl,yl); caxis([0,300])
% $$$ % $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ Z = [dsa(w==1&aper,'spine_lower','pelvis_root',2),...
% $$$      dsa(w==1&aper,'spine_lower','spine_middle',2),...
% $$$      dsa(w==1&aper,'spine_lower','head_front',3),...
% $$$      dsx(w==1&aper,'spine_lower',3),...
% $$$      abs(circ_dist(dsa(w==1&aper,'spine_lower','spine_upper',1),dsa(w==1&aper,'spine_middle','head_front',1)))];
% $$$ 
% $$$ 
% $$$ [State, hmm, decode] = gausshmm(nunity(Z),9);
% $$$ 
% $$$ 
% $$$ hmm =      hmmtrain(hmm.data.Xtrain,hmm.data.T,hmm);
% $$$ [decode] = hmmdecode(hmm.data.Xtrain,hmm.data.T,hmm);
% $$$ State = decode(1).q_star;
% $$$ 
% $$$ 
% $$$ ind = w==1&aper;
% $$$ Stateall = zeros(ufet.size(1),1);
% $$$ Stateall(ind.data) = State;
% $$$ figure
% $$$ plot(Stateall*10+1,'c')
% $$$ Lines(Trial.stc{'r',ufet.sampleRate}.data(:),[],'r');
% $$$ Lines(Trial.stc{'w',ufet.sampleRate}.data(:),[],'k');
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for i = 1:9,
% $$$     ind = w==1&aper&Stateall==i;
% $$$     %ind = w==1&aper;
% $$$     %x = dsa(w==1&aper&Stateall==i,'spine_lower','pelvis_root',2);   xl = .2:.01:1.6;
% $$$     %x = dsa(ind,'spine_lower','spine_middle',2);          xl = 0:.01:1.2;
% $$$     x = abs(circ_dist(dsa(ind,'spine_lower','spine_upper',1),dsa(ind,'spine_middle','head_front',1))); xl = 0:.02:pi;
% $$$     %x = ufet(w==i&aper);                                  xl = 18:.5:65;
% $$$     %x = dsx(w==i&aper,'spine_lower',3);                   xl = 1:1:70;
% $$$     y = dsx(ind,'spine_lower',3);                          yl = 1:1:70;
% $$$     %y = dsx(w==1&aper&Stateall==i,'pelvis_root',3);       yl = 15:1:110;
% $$$     %y = vel(w==1&aper&Stateall==i,'head_front');          yl = -1.5:.05:2;
% $$$     subplotfit(i,9);
% $$$     hist2([x,y],xl,yl); caxis([0,200])
% $$$ end
