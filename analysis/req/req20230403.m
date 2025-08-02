% req20230403
%    Tags: ccg synaptic transmition
%    Status: Active
%    Type: Utility
%    Author: Justin Graboski
%    Final_Forms: N/A
%    Project: General
%    Description: investigate synaptic transmition gain as function of theta phase


MjgER2016_load_data();

Trial = Trials{18}

unitInts = select_units(Trial,'int');
unitPyrs = select_units(Trial,'pyr');

ints = Trial.load('spk', Trial.lfp.sampleRate, 'theta-groom-sit-rear', unitInts);
pyrs = Trial.load('spk', Trial.lfp.sampleRate, 'theta-groom-sit-rear', unitPyrs);


phz = load_theta_phase(Trial,Trial.lfp.sampleRate);
phzBin = linspace(0,2*pi,11);

binSize = 4;
halfBins = 60;
mccg = zeros([121,10,numel(unitInts),numel(unitPyrs)]);
for int = unitInts
    for pyr = unitPyrs
        resI = ints(int);        
        resP = pyrs(pyr);
        for p = 1:numel(phzBin)-1
            %mresI = resI(phz(resI) >= phzBin(p) & phz(resI) < phzBin(p+1));
            mresI = resI;
            mresP = resP(phz(resP) >= phzBin(p) & phz(resP) < phzBin(p+1));            
            %mresP = resP;
            [tccg,timebins] =                                                    ...
                CCG([mresI;mresP],                                                 ...
                    [ones(size(mresI));2*ones(size(mresP))],                       ...
                    binSize,halfBins,                                            ...
                    Trial.lfp.sampleRate,                                        ...
                 [1,2],                                                       ...
                    'count');
            mccg(:,p,int==unitInts,pyr==unitPyrs) = tccg(:,1,2);
        end
    end
end


ii=9;
figure,
for p = 1:numel(unitPyrs)
    clf();
    subplot(311);
    imagesc(timebins,1:10,mccg(:,:,ii,p)');
    colorbar();
    subplot(312);
    imagesc(timebins,1:10,bsxfun(@rdivide,mccg(:,:,ii,p)',sum(mccg(:,:,ii,p)',2)));    
    colorbar();
    subplot(313);
    imagesc(timebins,1:10,bsxfun(@rdivide,mccg(:,:,ii,p)',sum(mccg(:,:,ii,p)',1)));    
    title(num2str(p));
    colorbar();
    waitforbuttonpress();
end
