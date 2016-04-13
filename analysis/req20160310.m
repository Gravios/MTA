function req20160310(Trial)
% PA - Heiarchical Segmentation of Behaviors 
% 1. Preproc features and Order selection
% 2. tsne dimensionallity reduction
% 3. Train neural networks on the target state vs all others
% 4. Accumulate stats of network output

Trial = 'jg05-20120317.cof.all';
Trial = MTATrial.validate(Trial);
%cd(Trial.path.data);
local = false;
overwrite = true;

% 1. Preproc features
file_preproc = fullfile(Trial.path.data,'analysis','req20160310_1_preproc.mat');
if (~exist(file_preproc,'file')&&~local)||overwrite,
    jid = popen(['MatSubmitLRZ --config lrzc_hugemem.conf'...
                 ' -y ' Trial.path.data ' -l ' Trial.name ' req20160310_1_preproc']);
    r1jid = [' -d afterok:' char(jid.readLine)];
else
    req20160310_1_preproc(Trial);
end

%2. t-SNE
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_2_tsne',num2str(s),'.mat');
    if (~exist(file_preproc,'file')&&~local)||overwrite,
        popen(['MatSubmitLRZ --config lrzc_hugemem.conf '...
               ' -y ' Trial.path.data r1jid ' -l ' Trial.name ...
               ' req20160310_2_tsne ',num2str(s)]);
    else
        req20160310_2_tsne(Trial,s);
    end
end

% 3. Train neural networks on the target state vs all others
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_3_trainNN',num2str(s),'.mat');
    if (~exist(file_preproc,'file')&&~local)||overwrite,
        jid = popen(['MatSubmitLRZ --config lrzc_mpp1.conf ' ...
                     ' -y ' Trial.path.data r1jid ' -l ' Trial.name ...
                     ' req20160310_3_trainNN ',num2str(s)]);
        r3jid = [' -d afterok:' char(jid.readLine)];
    else
        req20160310_3_trainNN(Trial,s);
    end
end

% 4. Accumulate stats of network output
for s = 1:5,
    file_preproc = fullfile(Trial.spath,'req20160310_4_accumStats',num2str(s),'.mat');
    if (~exist(file_preproc,'file')&&~local)||overwrite,
        popen(['MatSubmitLRZ --config lrzc_mpp1.conf ' ...
               ' -y ' Trial.path.data r3jid ' -l ' Trial.name ...
               ' req20160310_4_accumStats ',num2str(s)]);
    else
        req20160310_4_accumStats(Trial,s);
    end
end

%req20160310_5_genfigs
