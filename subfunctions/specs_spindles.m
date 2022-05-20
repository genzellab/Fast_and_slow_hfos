function [x,y,z,w,h,q,l,p]=specs_spindles(si,timeasleep,fn,print_hist)

    if ~isempty(si)

        %Instantaneous frequency.
        x=cellfun(@(equis) mean(instfreq(equis,fn)) ,si,'UniformOutput',false);
        x=cell2mat(x);
        if print_hist==1
            subplot(3,2,1)
%             xlim([100 250])
            histogram(x,20); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count')
%             ylim([0 35])
%             yticks([0:5:35])

%             histogram(x,[100:10:250],'Normalization','probability'); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.5])
%             yticks([0:0.1:0.5])
        end
        x=median(x);
        %Average frequency
        y=cellfun(@(equis) (meanfreq(equis,fn)) ,si,'UniformOutput',false);
        y=cell2mat(y);
 
        if print_hist==1
            subplot(3,2,2)
%             xlim([100 250])
            histogram(y,20); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
%             ylim([0 50])
 xlim([8 20])
 xticks([8:2:20])
%             yticks([0:10:50])
%             histogram(y,[100:10:250],'Normalization','probability'); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.5]) 
%             yticks([0:0.1:0.5])            
        end
        y=median(y);

        %Amplitude
        z=cellfun(@(equis) max(abs(hilbert(equis))) ,si,'UniformOutput',false);
        z=cell2mat(z);
        if print_hist==1
            subplot(3,2,3)
%             xlim([0 30])
            histogram(z,20); title('Amplitude');xlabel('\muV');ylabel('Count')         
%             ylim([0 60])
%             yticks([0:10:60])
%             histogram(z,[0:1:30],'Normalization','probability'); title('Amplitude');xlabel('\muV');ylabel('Probability')
%             ylim([0 0.5])
%             yticks([0:0.1:0.5])            
        end
        z=median(z);
        
        %Area under the curve
        l=cell2mat(cellfun(@(equis) trapz((1:length(equis))./1000,abs(equis)),si,'UniformOutput',false));
        if print_hist==1
            subplot(3,2,4)
%             xlim([0 0.5])
            histogram(l,20); title('Area under the curve');xlabel('AUC');ylabel('Count') 
%             ylim([0 50])
%             yticks([0:10:50])
%             histogram(l,[0:0.025:0.5],'Normalization','probability'); title('Area under the curve');xlabel('AUC');ylabel('Probability')
%             ylim([0 0.5])
%             yticks([0:0.1:0.5])            
        end
        l=median(l);

        %Count
        w=length(si);

        %Rate
        h=w/(timeasleep*(60));

        %Duration
        q=(cellfun('length',si)/1000);
        if print_hist==1
            subplot(3,2,5)
%             xlim([0 0.1]*1000)
%             xticks([0 0.025 0.05 0.075 0.1]*1000)
            histogram(q*1000,20); title('Duration');xlabel('MIliseconds');ylabel('Count')    
%             ylim([0 40])
%             yticks([0:10:40])
%             histogram(q*1000,[0:0.005:0.1].*1000,'Normalization','probability'); title('Duration');xlabel('MIliseconds');ylabel('Probability')
%             ylim([0 0.5])
%             yticks([0:0.1:0.5])
        end
        q=median(q);
        
        %Peak-to-peak distance
        p=cellfun(@peak2peak,si);
        if print_hist==1
            subplot(3,2,6)
%             xlim([0 45])
            histogram(p,20); title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Count');   
%             ylim([0 70])
%             yticks([0:10:70])
%             histogram(p,[0:2:45],'Normalization','probability'); title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Probability');
%             ylim([0 0.5]) 
%             yticks([0:0.1:0.5])            
        end
        p=median(p);


    else

        x=NaN;
        y=NaN;
        z=NaN;
        w=0;
        h=NaN;
        q=NaN;
        l=NaN;
        p=NaN;

    end

end