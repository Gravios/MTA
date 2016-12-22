function figId = gen_figure_id
figId = round((now-floor(now*1e-3)*1e3)*1e5);
