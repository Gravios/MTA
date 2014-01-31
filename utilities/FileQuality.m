% FileQuality(FileBase, ElecNo, BurstTimeWin, Fets2Use, Refrac,Verbose)
%
% a wrapper function that runs ClusterQuality for every cluster in a file.
%
% BurstTimeWin defaults to 120 (6 ms at 20 kHz)
% Fets2Use defaults to 1:12
% results are printed to the console
%
% optional output Out is a 4-column array with 1 row per cell
% (starting from 2) with columns CluNo, eDist, bRat, fraction of ISIs < Refrac
% (default 40, i.e. 2ms @ 20kHz);

function fq = FileQuality(FileBase, ElecNo, BurstTimeWin, Fets2Use, Refrac,Verbose)

Par = LoadPar([FileBase '.xml']);
[Fet nFet] = LoadFet([FileBase '.fet.' num2str(ElecNo)]);
Clu = LoadClu([FileBase '.clu.' num2str(ElecNo)]);
Res = load([FileBase '.res.' num2str(ElecNo)]);

if nargin<3, BurstTimeWin = 120*Par.SampleRate/1000*20; end;
if nargin<4, 
    %Par1 = LoadPar1([FileBase '.par.' num2str(ElecNo)]);
    %Fets2Use = 1:(Par1.nSelectedChannels*Par1.nPCs);
    Fets2Use = [1:nFet-1];
end
if nargin<5, Refrac = Par.SampleRate*2/1000; end;
if nargin<6 Verbose=1; end

Fet = Fet(:,Fets2Use);



for CluNo = 2:max(Clu)
    if sum(Clu==CluNo)>1
        [eDist(CluNo-1), bRat(CluNo-1)] = ...
            ClusterQuality(Fet, find(Clu==CluNo), Res, BurstTimeWin,Verbose);
%%%%%%%%%%debug
% $$$ if eDist(CluNo-1)==0 
% $$$     keyboard
% $$$ end
%%%%%%%%%%%%%%%
        MyISI = diff(Res(find(Clu==CluNo)));
        RefracViol(CluNo-1) = sum(MyISI<Refrac)/length(MyISI);
        [MaxFet(CluNo-1) mi] = max(mean(Fet(Clu==CluNo,1:3:end-1)));
        MaxFet(CluNo-1) = (MaxFet(CluNo-1) - mean(Fet(:,(mi-1)*3+1)))./std(Fet(:,(mi-1)*3+1));
    else
        RefracViol(CluNo-1) = 0;
        eDist(CluNo-1)=0;
        bRat(CluNo-1)=0;
        MaxFet(CluNo-1)=0;
    end
end


CluIds = unique(Clu);
CluInds = CluIds(~ismember(CluIds,[0,1]));


if nargout>=1
    fq = [CluInds' ; eDist(CluInds-1) ; bRat(CluInds-1); RefracViol(CluInds-1); MaxFet(CluInds-1)]';
    save([FileBase '.FileQuality.' num2str(ElecNo) '.mat'],'fq');
else
    fprintf('Cell %d: eDist %f bRat %f RefViol %f maxFet %f\n', ...
        [CluInds' ; eDist(CluInds-1) ; bRat(CluInds-1); RefracViol(CluInds-1); MaxFet(CluInds-1)]);
end