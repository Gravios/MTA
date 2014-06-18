function gausswindow = gtwin(twin,sampleRate)
    gausswindow = round(twin*sampleRate)+(mod(round(twin*sampleRate),2)==0);    
    gausswindow = gausswin(gausswindow)./sum(gausswin(gausswindow));
end