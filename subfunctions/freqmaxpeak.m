function[freqmax] =freqmaxpeak(x,fn)
[power,freq]=periodogram(x,[],[],fn);

[~,i2]=max(((sqrt(power))).*freq);
freqmax=freq(i2);

end