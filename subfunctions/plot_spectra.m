function plot_spectra(P,Q,labelconditions2,label1,s,w,same_nr_types,N)

if same_nr_types==1
    n=N;
else
    n=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);    
end


    %Order ripples
    p=P.plusmaze.(label1{w}){s}; 
    q=Q.plusmaze.(label1{w}){s}; 
    % R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,:)))),q));
    R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,121-50:121+50)))),q));
    [~,r]=sort(R,'descend');
    p=p(r);
    q=q(r);
    p=p(1:n);
    q=q(1:n);

    %Order ripples
    p_nl=P.nl.(label1{w}){s}; 
    q_nl=Q.nl.(label1{w}){s}; 
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

    if w==3 %PAR-centered ripples.
         p=cellfun(@(equis1) flip(equis1),p,'UniformOutput',false);
         q=cellfun(@(equis1) flip(equis1),q,'UniformOutput',false);
         p_nl=cellfun(@(equis1) flip(equis1),p_nl,'UniformOutput',false);
         q_nl=cellfun(@(equis1) flip(equis1),q_nl,'UniformOutput',false);
    end
    
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
    g=title([label1{j} ' Baseline']);
    g.FontSize=12;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    ylim([100 250])

    subplot(3,3,3*j-1)
    ft_singleplotTFR(cfg, freq2); 
    g=title(strcat([label1{j} ' Plusmaze' ]));
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

    g=title(strcat(labelconditions2{1},' vs Baseline'));
    g.FontSize=12;
    ylim([100 250])
    

    end

end
