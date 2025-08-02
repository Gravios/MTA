function req20160310(Trial)
% PA - Heiarchical Segmentation of Behaviors 
% 1. Preproc features and Order selection
% 2. tsne dimensionallity reduction
% 3. Train neural networks on the target state vs all others
% 4. Accumulate stats of network output

%Trial = 'jg05-20120317.cof.all';
Trial = 'Ed03-20140624.cof.all';
Trial = MTATrial.validate(Trial);
%cd(Trial.path.project);
local = true;
overwrite = true;
sampleRate = 12;
states = {'walk','rear','turn','pause','groom','sit'};

% PREPROCCESS (1) features and save to relevant files for subsequent analysis
file_preproc = fullfile(Trial.path.project,'analysis','req20160310_1_preproc.mat');
if ~exist(file_preproc,'file')||overwrite
    if ~local,
        jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf'...
                                  ' -l ' Trial.name  ...
                                  ' req20160310_1_preproc']);
        r1jid = [' -d afterok:' char(jid.readLine)];
    else
        req20160310_1_preproc(Trial,sampleRate,states,true,refTrial,refStcMode);
    end
end

% TSNE (2) t-SNE
r1jid='';
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_2_tsne',num2str(s),'.mat');
    if ~exist(file_preproc,'file')||overwrite,
        if ~local,
            system(['MatSubmitLRZ --config lrzc_hugemem.conf '...
                    r1jid ' -l ' Trial.name ...
                    ' req20160310_2_tsne ',num2str(s)]);
        else
            req20160310_2_tsne(Trial,s);
        end
    end
end


% TRAIN (3) neural networks on the target state vs all others
r1jid='';
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_3_trainNN',num2str(s),'.mat');
    if ~exist(file_preproc,'file')||overwrite,
        if ~local,
            jid = popen(['MatSubmitLRZ --config lrzc_serial.conf ' ...
                         r1jid ...
                         ' -l ' Trial.name ...
                         ' req20160310_3_trainNN '...
                         num2str(s)]);
            r3jid = [' -d afterok:' char(jid.readLine)];
        else
            req20160310_3_trainNN(Trial,s);
        end
    end
end


% ACCUMULATE (4) stats of network output
r3jid='';
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_4_accumStats',num2str(s),'.mat');
    if (~exist(file_preproc,'file')||overwrite)&&~local,
        system(['MatSubmitLRZ --config lrzc_serial.conf ' ...
                              r3jid ' -l ' Trial.name ...
                              ' req20160310_4_accumStats ',num2str(s)]);
    elseif (~exist(file_preproc,'file')||overwrite)
        req20160310_4_accumStats(Trial,s);
    end
end


% GENERATE (5) NN figures
%    a. sort feature order to maximize accuracy gain from incremental
%       addition of features to neural network model. 
%    b. generate figures for suplementary plots
%
file_preproc = fullfile(Trial.spath,'req20160310_5_genfigs.mat');
if ~exist(file_preproc,'file')||overwrite,
    req20160310_5_genfigs(Trial)
end


% TRAIN (6) NN with optimized feature order
%    train nn's with incremental addition of features selected
%    during req20160310_5_genfigs

r5jid='';
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_6_trainOptNN',num2str(s),'.mat');
    if ~exist(file_preproc,'file')||overwrite,
        if ~local,
            jid = popen(['MatSubmitLRZ --config lrzc_serial.conf ' ...
                         r5jid ...
                         ' -l ' Trial.name ...
                         ' req20160310_6_trainOptNN '...
                         num2str(s)]);
            r6jid = [' -d afterok:' char(jid.readLine)];
        else
            req20160310_6_trainOptNN(Trial,s);
        end
    end
end


% ACCUMULATE (7) stats of optimized network output
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_7_accumOptStats',num2str(s),'.mat');
    if (~exist(file_preproc,'file')||overwrite)
        if ~local,
            jid = popen(['MatSubmitLRZ --config lrzc_serial.conf ' ...
                         r6jid ' -l ' Trial.name ...
                         ' req20160310_7_accumOptStats ',num2str(s)]);
            r6jid = [' -d afterok:' char(jid.readLine)];    
        else (~exist(file_preproc,'file')||overwrite)
            req20160310_7_accumOptStats(Trial,s);
            
        end
    end
end


% GENERATE (8) Optimized NN figures
file_preproc = fullfile(Trial.spath,'req20160310_8_genfigs.mat');
if ~exist(file_preproc,'file')||overwrite,
    req20160310_5_genfigs(Trial);
end


