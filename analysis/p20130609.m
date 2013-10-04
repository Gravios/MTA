function p20130609()
%batch_trial_list,batch_tag,states,trialname)

trialname = 'all'
batch_trial_list = '/data/homes/gravio/data/analysis/batch_trial_list';
batch_tag = 'rear_pfc';
states = {'walk','rear','theta'};

fid = fopen(batch_trial_list);

tline = fgets(fid);
while ischar(tline)
    for i = 1:length(states),
        consolidate_batch(tline(1:end-1),trialname,[batch_tag '_' states{i}],1,0)
    end
    tline = fgets(fid);
end


% $$$ system(['scp -v knajg01@hpc-bw.uni-tuebingen.de:/home-link/knajg01/data/' ...
% $$$     'analysis/ ' Trial.name '/' Trial.filebase '.ccg.' batch_tag '* ' Trial.spath.analysis])

