function quick_qda(Trial,varargin)
[New_Stc_Mode,QDA_model,debug,train] = DefaultArgs(varargin,{'qda_filtf1p5','MTA_mknsrw_filtf1p5_QDA_model.mat',false,false});


if train,
    Trial.stc.updateMode('manual_mknsrw');
    Trial.stc.load;
    bhv_qda(Trial,[],true,[],QDA_model);
end

Stc = Trial.stc.copy;
Stc.updateMode(New_Stc_Mode);
Stc.states={};
Stc = bhv_qda(Trial,Stc,false,[],QDA_model);
if debug,keyboard,end
Stc.save(1);
