function pfbrowser(Session,varargin)
% function pf3viewer(Session,varargin)
[sgrid,display_ccg,mazeName,trialName]= DefaultArgs(varargin,{struct('x','trial','y','state','state',{{'rear','rear.theta','walk','walk.theta'}},'trial',{{'crt1','alt1','crt2','alt2'}}),1,1,'cof',''});

if ~isa(Session,'MTASession'),
    Session = MTASession(Session,mazeName);
end


statefield = {};
stateCCG = {};
numXticks = length(getfield(sgrid,sgrid.x));
numYticks = length(getfield(sgrid,sgrid.y));
numT = [];
numS = [];
if strcmp(sgrid.x,'trial'), 
    numT = numXticks;
    numS = numYticks;
else 
    numT = numYticks;
    numS = numXticks;
end

for i = 1:numXticks,
    for t = 1:numYticks,
        if strcmp(sgrid.x,'trial'), 
            tn = i;
            sn = t;
        else
            tn = t;
            sn = i;
        end
        
        Trial = MTATrial(Session,sgrid.trial{tn});
        if exist([Trial.spath.analysis Trial.filebase '.pf3.' Trial.trackingMarker '.' sgrid.state{sn} '.' num2str(numZslices) '.mat'],'file'),
            statefield{sn,tn} = load([Trial.spath.analysis Trial.filebase '.pf3.' Trial.trackingMarker '.' sgrid.state{sn} '.' num2str(numZslices) '.mat']);
            if display_ccg,
                stateCCG{sn,tn} = load([Trial.spath.analysis Trial.filebase '.ccg.' sgrid.state{sn} '.mat']);
            end
        end
    end
end

if exist([Session.spath.nlx Session.name '.NeuronQuality.mat'],'file'),
    load([Session.spath.nlx Session.name '.NeuronQuality.mat']);
    [~, si] = sort(nq.eDist,'descend');
    unit=si(1);
end




while 1
    figure(100)
    set(gcf,'Name',num2str(statefield{1,1}.Map(unit,1)'))
    clf
    if display_ccg
        ny = (numZslices-1+numYticks)*2+numYticks;
    else
        ny = (numZslices-1+numYticks)*2;
    end
    nx = length(unit)*numXticks*2;
    if nq.FirRate(unit)>0.2,
        units = unit;

        for i = 0:length(getfield(sgrid,sgrid.x))-1
            for t = 0:length(getfield(sgrid,sgrid.y))-1
                if strcmp(sgrid.x,'trial'), 
                    tn = i+1;
                    sn = t+1;
                else
                    tn = t+1;
                    sn = i+1;
                end
                for u= 1:1,
                    unit = units(1);
                    for slice = 1:numZslices,

                        subplot2(ny,nx,[t*ny/numYticks+1,t*ny/numYticks+2],[i*2+u*2-1,i*2+u*2]);
                        imagesc(statefield{sn,tn}.bin1{unit,slice},statefield{sn,tn}.bin2{unit,slice},statefield{sn,tn}.rateMap{unit,slice}')
                        colorbar
                        colormap('default')
                        cc = colormap;
                        cc(1,:) = [0 0 0];
                        colormap(cc)
                        imagesc(statefield{sn,tn}.bin1{unit,slice},statefield{sn,tn}.bin2{unit,slice},statefield{sn,tn}.rateMap{unit,slice}')
                        text(-450,200,sprintf('%2.1f',max(statefield{sn,tn}.rateMap{unit,slice}(:))),'Color','w','FontWeight','bold','FontSize',10);
                        axis xy
                    end

                    if display_ccg,
                        subplot2(ny,nx,t*3+3,i*2+u*2-1);
                        bar(stateCCG{sn,tn}.tbin, stateCCG{sn,tn}.sgccg(:,1,unit));axis tight; grid on
                        xlim([-2500 2500]);
                        subplot2(ny,nx,t*3+3,i*2+u*2);
                        bar(stateCCG{sn,tn}.tbin, stateCCG{sn,tn}.sgccg(:,2,unit));axis tight; grid on
                        xlim([-2500 2500]);
                    end

                    ForAllSubplots('set(gca,''FontSize'',10)');
                    if tn==1&sn==1,
                        fprintf('Unit# %d, Electrode %d, Cluster %d\n',unit,statefield{1,1}.Map(unit,2),statefield{1,1}.Map(unit,3));
                        fprintf('Quality: %1.2f, Refrac: %1.3f, SpkWidthR %1.1f msec, FirRate %2.1f Hz\n',nq.eDist(unit),nq.Refrac(unit),nq.SpkWidthR(unit),nq.FirRate(unit));
                    end
                end
            end
        end

        if length(units)>1 unit=units(end);end    
        wb =  waitforbuttonpress
        whatkey = get(gcf,'CurrentCharacter');
    end
    if ~wb, whatkey = 'n';end

    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
      case 30
        sii = find(si==unit);
        unit = si(max(1,sii-1));
      case 31
        sii = find(si==unit);
        unit = si(min(size(statefield{1,1}.Map,1),sii+1));
      case double('s') %save
        print('-dpng', '-r300', ['/tmp/plf-' num2str(units) '.png']);
      case double('b')
        keyboard
    end

    if unit<1  unit=1; end
    if unit>size(statefield{1,1}.Map,1) unit=1; end
end

    
