rat = load_patch_model('rat');
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false);
subject = update_subject_patch(subject,'body',[],false);
