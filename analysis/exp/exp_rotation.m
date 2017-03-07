s = MTASession('jg05-20120315');
Trial = QuickTrialSetup(s,'ctl1',[10,0],[2:7]);
Trial = QuickTrialSetup(s,'alt1',[10,0],[1,3:7]);
Trial = QuickTrialSetup(s,'ctl2',[10,0],[1:2,4:7]);
Trial = QuickTrialSetup(s,'alt2',[10,0],[1:3,5:7]);
Trial = QuickTrialSetup(s,'ctl3',[10,0],[1:4]);

Trial = MTATrial('jg05-20120315');
stc = Trial.stc.copy;

stc.updateMode('hand_labeled');
stc.load(Trial);

stc.updateMode('auto_wbhr');
stc.load(Trial);

Trial = MTATrial('jg05-20120315','ctl1');
Trial = MTATrial('jg05-20120315','alt1');
Trial = MTATrial('jg05-20120315','ctl2');
Trial = MTATrial('jg05-20120315','alt2');
Trial = MTATrial('jg05-20120315','ctl3');

trials = {'ctl1','alt1','ctl2','alt2','ctl3'};
for s = 1:numel(trials)
Trial = MTATrial('jg05-20120315',trials{s});    
stc.updateMode('auto_wbhr');
stc.load(Trial);
Trial.stc = stc;

display = true;


%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'theta','rear','walk','hang','lang','hswalk','lswalk'};
%states = {'hswalk&theta','lswalk&theta'};
%states = {'theta'};
nsts = size(states,2);
mode = 'std';
niter = 10000;


overwrite = true;
units = select_units(Trial,18,'pyr');

switch mode
  case 'std'
    for i = 1:nsts,
        MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
    end

  case 'shuff'
    for i = 1:nsts,
        MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
    end
end

end

if display,
    pfs ={};
    for s = 1:numel(trials)
Trial = MTATrial('jg05-20120315',trials{s});    
stc.updateMode('auto_wbhr');
stc.load(Trial);
Trial.stc = stc;


    switch mode
      case 'std'
        for i = 1:nsts,
            %i = 1;

            pfs{i,s} = MTAAknnpfs(Trial,units,states{i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
        
      case 'shuff'     
        for i = 1:nsts,
            pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
    
    end
end
    
    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for i = 1:nsts,
        for s = 1:numel(trials),
            subplot2(nsts,numel(trials),i,s);cla
            pfs{i,s}.plot(unit,[],true);
            title([pfs{i,s}.parameters.states,': ',num2str(unit)]);

        end

        end
          unit = figure_controls(hfig,unit,units);
    end
    

end

