function [n]=plot_spectra_ripples_vs_hfos(P,Q,SP,SQ,labelconditions2,label1,same_nr_types,N,str_P,str_SP)

if same_nr_types==1
    n=N;
else
    n=min([length([([ P.('plusmaze').('HPC'){2}])... 
           ([ P.('nl').('HPC'){2}])...
           ([ P.('for').('HPC'){2}])...
           ([ P.('novelty').('HPC'){2}])...
         ])...
        
        length([([ SP.('plusmaze').('PAR'){1} SP.('plusmaze').('PAR'){2}])... 
           ([SP.('nl').('PAR'){1} SP.('nl').('PAR'){2}])...
           ([SP.('for').('PAR'){1} SP.('for').('PAR'){2}])...
           ([SP.('novelty').('PAR'){1} SP.('novelty').('PAR'){2}])...
         ]) 
         
         
         ]);    
end


    %Order HPC ripples
    p=[P.plusmaze.('HPC'){2} P.nl.('HPC'){2} P.for.('HPC'){2} P.novelty.('HPC'){2}];
    q=[Q.plusmaze.('HPC'){2} Q.nl.('HPC'){2} Q.for.('HPC'){2} Q.novelty.('HPC'){2}]; 
    % R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,:)))),q));
    R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,121-50:121+50)))),q));
    [~,r]=sort(R,'descend');
    p=p(r);
    q=q(r);
    p=p(1:n);
    q=q(1:n);

    %Order (slow or fast) HFOs ripples
    p_nl=[SP.('plusmaze').('PAR'){1} SP.('plusmaze').('PAR'){2}...
          SP.('nl').('PAR'){1} SP.('nl').('PAR'){2} ...
          SP.('for').('PAR'){1} SP.('for').('PAR'){2} ...
          SP.('novelty').('PAR'){1} SP.('novelty').('PAR'){2}
        ]; 
    
    q_nl=[SQ.('plusmaze').('PAR'){1} SQ.('plusmaze').('PAR'){2}...
          SQ.('nl').('PAR'){1} SQ.('nl').('PAR'){2} ...
          SQ.('for').('PAR'){1} SQ.('for').('PAR'){2} ...
          SQ.('novelty').('PAR'){1} SQ.('novelty').('PAR'){2}
        ]; 
        
    % R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,:)))),q_nl));
    R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,121-50:121+50)))),q_nl));
    [~,r_nl]=sort(R,'descend');
    p_nl=p_nl(r_nl);
    q_nl=q_nl(r_nl);
    p_nl=p_nl(1:n);
    q_nl=q_nl(1:n);


    %Max 1000 ripples.
    if length(q)>1000
        q=q(1:1000);
        p=p(1:1000);
        q_nl=q_nl(1:1000);
        p_nl=p_nl(1:1000);
    end

%     if w==3 %PAR-centered ripples.
%          p=cellfun(@(equis1) flip(equis1),p,'UniformOutput',false);
%          q=cellfun(@(equis1) flip(equis1),q,'UniformOutput',false);
         p_nl=cellfun(@(equis1) flip(equis1),p_nl,'UniformOutput',false);
         q_nl=cellfun(@(equis1) flip(equis1),q_nl,'UniformOutput',false);
%     end

    
    ro=150;

    toy=[-.1:.001:.1];
    freq1=time_frequency(q_nl,create_timecell(ro,length(q_nl)),[100:1:300],[],toy);
    freq2=time_frequency(q,create_timecell(ro,length(q)),[100:1:300],[],toy);

allscreen()
    for j=1:3

    cfg              = [];
    cfg.channel      = freq1.label{j};
    [ zmin1, zmax1] = ft_getminmax(cfg, freq1);
    [zmin2, zmax2] = ft_getminmax(cfg, freq2);
    zlim=[min([zmin1 zmin2]) max([zmax1 zmax2])];

    cfg              = [];
    cfg.zlim=zlim;
    cfg.channel      = freq2.label{j};
    cfg.colormap=colormap(jet(256));

    subplot(3,3,3*j-2)
    ft_singleplotTFR(cfg, freq1); 
    g=title([label1{j} ' during ' str_SP]);
    g.FontSize=12;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    ylim([100 250])

    subplot(3,3,3*j-1)
    ft_singleplotTFR(cfg, freq2); 
    g=title(strcat([label1{j} ' during ' str_P ]));
    g.FontSize=12;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    ylim([100 250])


    % Pixel-based stats
    zmap=stats_high(freq1,freq2,j);
    subplot(3,3,3*j);

    colormap(jet(256))
    zmap(zmap == 0) = NaN;
    J=imagesc(freq1.time,freq1.freq,zmap)
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    set(gca,'xlim',xlim,'ydir','no')
    set(J,'AlphaData',~isnan(zmap))
    c=narrow_colorbar()
     c.YLim=[-max(abs(c.YLim)) max(abs(c.YLim))];
    caxis([-max(abs(c.YLim)) max(abs(c.YLim))])
    c=narrow_colorbar()

    g=title([label1{j} ' (' str_P,' vs ' str_SP ')']);
    g.FontSize=12;
    ylim([100 250])
    

    end

end
