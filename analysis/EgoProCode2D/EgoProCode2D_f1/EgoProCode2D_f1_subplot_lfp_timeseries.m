% EgoProCode2D_f1_subplot_lfp_timeseries

[pfig, sax] = setup_figure_( mfilename() );

timeIndex = 87000;
lfp = Trials{exampleLFP.trialIndex}.load('lfp',Trials{exampleUnit.trialIndex}.meta.channelGroup.theta);
lfp.resample(sampleRate);
lfp.data = nunity(lfp.data);
phz = load_theta_phase(Trials{exampleLFP.trialIndex},sampleRate);
mres = spk{exampleLFP.trialIndex}(exampleLFP.unitId); 
mres  = mres(WithinRanges(mres,timeIndex+[-200,200]));


plot(([-200:200])./sampleRate,lfp(timeIndex+[-200:200],1)+3,'k','LineWidth',1)
plot(([-200:200])./sampleRate,phz(timeIndex+[-200:200],1),'g','LineWidth',1)
scatter((mres-timeIndex)./sampleRate,lfp(mres,1)+3,4,'r','filled');
sax(end).Visible = 'off';

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

