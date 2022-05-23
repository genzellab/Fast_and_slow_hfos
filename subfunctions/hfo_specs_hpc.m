function [x,y,z,w,h,q,l,p,si_mixed,th]=hfo_specs_hpc(si,timeasleep,print_hist,fn)

    if ~isempty(si)

        %Instantaneous frequency.
        x=cellfun(@(equis) mean(instfreq(equis,fn)) ,si,'UniformOutput',false);
        x=cell2mat(x);
        if print_hist==1
            subplot(3,2,1)
            xlim([100 250])
            histogram(x,[100:2:250]); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count')
            xticks([100:50:250])
%             ylim([0 500])
%             yticks([0:100:500])

%             histogram(x,[100:2:250],'Normalization','probability'); title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])
        end
        x=median(x);
        %fi_cortex(k)=x;
        %Average frequency
        y=cellfun(@(equis) (meanfreq(equis,fn)) ,si,'UniformOutput',false);
        y=cell2mat(y);
%         th=gaussmix(y);
%         si_mixed.g1=si(y<=th);
%         si_mixed.i1=find(y<=th);
%         si_mixed.g2=si(y>th);
%         si_mixed.i2=find(y>th);
th=[];
si_mixed=[];
 
        if print_hist==1
            subplot(3,2,2)
            xlim([100 250])
            xticks([100:50:250])

            histogram(y,[100:2:250]); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
%             ylim([0 700])
%             
%             yticks([0:100:700])
%             histogram(y,[100:2:250],'Normalization','probability'); title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Probability')
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])            
        end
        y=median(y);
        %fa_cortex(k)=y;

        %Amplitude
        z=cellfun(@(equis) max(abs(hilbert(equis))) ,si,'UniformOutput',false);
        z=cell2mat(z);
        if print_hist==1
            subplot(3,2,3)
            histogram(z,[0:2:120]); title('Amplitude');xlabel('\muV');ylabel('Count')         
%             ylim([0 1000])
%             yticks([0:250:1000])
            xlim([0 120])
            xticks([0:30:120])

%             histogram(z,[0:2:120],'Normalization','probability'); title('Amplitude');xlabel('\muV');ylabel('Probability')
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])            
        end
        z=median(z);
        %amp_cortex(k)=z;
        
        %Area under the curve
        l=cell2mat(cellfun(@(equis) trapz((1:length(equis))./fn,abs(equis)),si,'UniformOutput',false));
        if print_hist==1
            subplot(3,2,4)
            histogram(l,[0:0.10:4]); title('Area under the curve');xlabel('AUC');ylabel('Count') 
%             ylim([0 1500])
%             yticks([0:250:1500])
            xlim([0 4])
            xticks([0:0.5:4])
%             histogram(l,[0:0.10:4],'Normalization','probability'); title('Area under the curve');xlabel('AUC');ylabel('Probability')
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])            
        end
        l=median(l);

        %Count
         w=length(si);
%w=[];

        %Rate
         h=w/(timeasleep*(60));
%h=[];

        %Duration
        q=(cellfun('length',si)/fn);
        if print_hist==1
            subplot(3,2,5)

            histogram(q*1000,[0:0.01:0.6].*1000); title('Duration');xlabel('MIliseconds');ylabel('Count')    
%             ylim([0 2000])
%             yticks([0:500:2000])
            xlim([0 0.4]*1000)
            xticks([0:0.1:0.4]*1000)
%             histogram(q*1000,[0:0.01:0.6].*1000,'Normalization','probability'); title('Duration');xlabel('MIliseconds');ylabel('Probability')
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])
        end
        q=median(q);
        
        %Peak-to-peak distance
        p=cellfun(@peak2peak,si);
        if print_hist==1
            subplot(3,2,6)
           
            histogram(p,[0:5:200]); title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Count');   
%             ylim([0 1500])
%             yticks([0:500:1500])
            xlim([0 200])
            xticks([0:50:200])
%             histogram(p,[0:5:200],'Normalization','probability'); title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Probability');
%             ylim([0 0.2])
%             yticks([0:0.05:0.2])            
        end
        p=median(p);


    else

        x=NaN;
        y=NaN;
        z=NaN;
        w=NaN;
        h=NaN;
        q=NaN;
        l=NaN;
        p=NaN;

    end

end