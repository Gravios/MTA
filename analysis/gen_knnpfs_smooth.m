%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
%MTAstartup;
%Trial = MTATrial('er06-20130614','fly');
%Session = MTASession('er06-20130614');
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');
display = true;

%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
%states = {'hswalk&theta','lswalk&theta'};
states = {'walk'};
nsts = size(states,2);
mode = 'shuff';
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

if display,
    pfs ={};
    switch mode
      case 'std'
        for i = 1:nsts,
            Pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
        
      case 'shuff'     
        for i = 1:nsts,
            pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
    
    end

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for i = 1:nsts,
            subplotfit(i,nsts);cla
            pfs{i}.plot(unit,[],true);
            title([pfs{i}.parameters.states,': ',num2str(unit)]);
        end
        unit = figure_controls(hfig,unit,units);
    end

end


hfig = figure;
unit = units(1);
while unit~=-1,
    for i = 1:nsts,
        subplot2(nsts,2,i,1);cla
        pfs{i}.plot(unit,[],true);
        subplot2(nsts,2,i,2);cla
        pfs{i}.plot(unit,[],true);
        title([pfs{i}.parameters.states,': ',num2str(unit)]);
    end
    unit = figure_controls(hfig,unit,units);
end

