%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
%MTAstartup;
%Trial = MTATrial('er06-20130614','fly');
%Session = MTASession('er06-20130614');
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');

Trial = MTATrial('Ed10-20140820','telcrtl1','rov');
Trial = MTATrial('Ed10-20140820','telshift1','rov');
Trial = MTATrial('Ed10-20140820','telcrtl2','rov');
Trial = MTATrial('Ed10-20140820','telshift2','rov');
display = true;

%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
%states = {'hswalk&theta','lswalk&theta'};
states = {'theta'};
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
        %for i = 1:nsts,
        i = 1;
        pfs{i} = MTAAknnpfs(Trial,units,states{1},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
        
      case 'shuff'     
        for i = 1:nsts,
            pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);            
        end
    
    end

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for i = 1:4,
            subplotfit(i,4);cla
            imagesc(reshape(pfs{i}.data.rateMap(:,unit),35,60)),colorbar
            %pfs{i}.plot(unit,'isCircular',false);
            title([pfs{i}.parameters.states,': ',num2str(unit)]);
        end
        unit = figure_controls(hfig,unit,units);
    end

end

