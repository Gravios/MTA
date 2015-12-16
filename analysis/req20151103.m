%% Neural network stuff
% 1. Train a neural network on the data of jg05-20120317. 
% 2. Behavior of other animals is then classified with the network. 
% 3. The Results are compared to the hand labeled data with a confussion matrix

states = {'walk','rear','turn','groom','pause','sit'};
newSampleRate = 15;

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');



x = fet_tsne(Trial,15,true);
%x = fet_all(Trial,15,true);
t = MTADxyz('data',stc2mat(Trial.stc,x,states),'sampleRate',x.sampleRate);
ind = Trial.stc{'a'};


net = patternnet(400);
view(net);
[net,tr] = train(net,x(ind,:)',~~t(ind,:)');


%Test other sessions 

Trial = MTATrial('er01-20110719');
Trial = MTATrial('Ed01-20140707');
Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev1');   
%Trial.stc.states{2}.data = Trial.stc.states{8}.data; %for er01
xyz   = Trial.load('xyz');
xyz.resample(newSampleRate);

xo = fet_tsne(Trial,15,true,mx,sx);
%xo = fet_all(Trial,newSampleRate,true);
ind = Trial.stc{'a'};
[y] = net(xo(ind,:)')';


shl = MTADxyz('data',double(0<stc2mat(Trial.stc,xyz,states)),'sampleRate',xyz.sampleRate);
shl = shl(ind,:);
ysm = bsxfun(@eq,y,max(y,[],2));

cm = confmat(shl,ysm); % DEP: netlab                
precision = round(diag(cm)'./sum(cm),4).*100;
sensitivity = round(diag(cm)./sum(cm,2),4).*100;
cm = round(cm./xyz.sampleRate,2);

figure(1),clf
plotmat([[cm,sensitivity];[precision,0]],'b','k',14);
h = gca;
h.XAxisLocation = 'top';
h.XTickLabel = cellfun(@horzcat,repmat({'     '},1,8),{states{:},'sensitivity',''},'uniformoutput',false)
%h.XTickLabel = {states{:},'sensitivity'};
tstates = {states{:},'precision',''};
%padding = repmat(' ',1,max(cellfun(@length,tstates)));
% $$$ tstates = cellfun(@vertcat,repmat({padding},1,numel(tstates)),cellfun(@horzcat,tstates,cellfun(@repmat,repmat({' '},1,numel(tstates)),repmat({1},1,numel(tstates)),mat2cell(abs(cellfun(@length,tstates)-max( cellfun(@length,tstates))),1,ones(1,numel(tstates))),'uniformoutput',false),'uniformoutput',false),'uniformoutput',false)

tstates = cellfun(@horzcat,tstates,cellfun(@repmat,repmat({' '},1,numel(tstates)),repmat({1},1,numel(tstates)),mat2cell(abs(cellfun(@length,tstates)-max( cellfun(@length,tstates))),1,ones(1,numel(tstates))),'uniformoutput',false),'uniformoutput',false);
h.YTickLabel = fliplr(tstates);
suptitle({['neural network vs hand labels for ',Trial.name],'nn trained jg05-20120317'})



Trial.fet = MTADfet(Trial.spath,...
                      [],...
                      [],...
                      [],...
                      Trial.sync.copy,...
                      Trial.sync.data(1),...
                      []);                  


Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev1'); % for er01: Trial.stc.states{2}.data = Trial.stc.states{8}.data;
Trial.stc.mode = 'hand_labeled_rev1';
Trial.stc.filename = 'Ed03-20140624.cof.all.stc.hand_labeled_rev1.mat';
Trial.stc.states{1}.key = 'p';
Trial.stc.states{3}.key = 'r';
Trial.stc.states{4}.key = 'n';
Trial.stc.states{5}.key = 'k';
Trial.stc.states{6}.key = 's';
Trial.stc.states{7}.key = 'm';
Trial.stc.save(1);

