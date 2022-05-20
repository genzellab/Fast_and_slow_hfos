
function [mdam,mdam2,mdam3,mdam4]=small_window(freq2,w,win_size)
%Compute mean power value of window for different frequency bands

dam=((squeeze(mean(squeeze(freq2.powspctrm(:,w,:,1+win_size:end-win_size)),1)))); %Average all events.
mdam=mean(dam(:)); %Mean value 

freqs=freq2.freq;

%100 to 150 Hz
n1=sum(freqs<=150);
dam=((squeeze(mean(squeeze(freq2.powspctrm(:,w,1:n1,1+win_size:end-win_size)),1)))); %Average all events.
mdam2=mean(dam(:)); %Mean value 

%151 to 200 Hz
n2=sum(freqs<=200);
dam=((squeeze(mean(squeeze(freq2.powspctrm(:,w,n1+1:n2,1+win_size:end-win_size)),1)))); %Average all events.
mdam3=mean(dam(:)); %Mean value 

%201 to 250 Hz
n3=sum(freqs<=250);
dam=((squeeze(mean(squeeze(freq2.powspctrm(:,w,n2+1:n3,1+win_size:end-win_size)),1)))); %Average all events.
mdam4=mean(dam(:)); %Mean value 

end

