function sniff_ccg(Session)

Session = MTASession('jg05-20120310');
Session = MTATrial(Session,'all');

Session = Session.load_ang();

dist = Session.ang(:,4,5,3);

dist(isnan(dist)) = mean(dist(~isnan(dist)));

fdist = Filter0(gausswin(5)./sum(gausswin(5)),dist);

wdist = WhitenSignal(fdist);


nffts = 2^7;
windows = 2^6;
frange = [1 20];

[y,f,t] = mtchglong(wdist,nffts,Session.xyzSampleRate,windows, ...
                    windows*.875,[],[],[],frange);

figure,imagesc(t+diff(t(1:2))/2,f,log10(y')),axis xy

sp_out = mean(y(:,f<6|f>14),2);
sp_in = mean(y(:,f>6&f<14),2);

spiors =sp_in-sp_out;
spior = sp_in./sp_out;


wdist = WhitenSignal(fdist);


sfs = 1/diff(t(1:2));

nffts = 2^6;
windows = 2^5;
frange = [1 20];

wspior = WhitenSignal(spior);

[ys,fs,ts] = mtchglong(wspior,nffts,sfs,windows,windows*.875,[],[],[],frange);

figure,imagesc(ts+diff(ts(1:2))/2,fs,log10(ys')),axis xy


sps_out = mean(ys(:,fs<4|fs>6),2);
sps_in = mean(ys(:,fs>4&fs<6),2);

figure,plot(ts,Filter0(gausswin(7)./sum(gausswin(7)),sps_in./sps_out))


fsps = Filter0(gausswin(7)./sum(gausswin(7)),sps_in./sps_out);

sres = round((ts(LocalMinima(-fsps,4,-2))+diff(ts(1:2))/2)*Session.xyzSampleRate);

sres = round((t(LocalMinima(-spiors,4,-0.002))+diff(t(1:2))/2)*Session.xyzSampleRate);

%% Batch Information
time_stamp = '20121205T1714';

%name - string: unique descriptive name of the ccg
name = 'rand_shit_sniff';

%Description - string: short statement describing the ccg
Description = 'oh god I do not know'

Res = {sres}
Clu = {1};
CluTags = {'rand_shit'};
method = 'normal';
surrogate_sample_size = 0;
partitions = 1;
numIterations = 1;%1000;

%% Check if Backups Exist
if ~exist([ Session.path.root 'batch/' name '-' time_stamp '.ccg.mat' ],'file'),
    bccg_search = MTAccg(Session,name,Description,{},Clu,CluTags,method,partitions,[],1,surrogate_sample_size,[],numIterations);
    %% Save a batch seed used to consolidate the results of a parallel batch
    save([ Session.path.root 'batch/' name '-' time_stamp '.ccg.mat' ],'bccg_search','-v7.3');
    %% Save a script backup 
    oriScript = mfilename('fullpath');
    system(['cp ' oriScript '.m ' Session.path.root 'batch/scripts/' mfilename '-' time_stamp '.m']);
end

tic
%% Run CCG
bccg = MTAccg(Session,name,Description,Res,Clu,CluTags,method,partitions,{},1,[],[],numIterations);
toc



