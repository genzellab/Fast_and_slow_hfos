function plot_spec_traces(P,Q,labelconditions2,label1,s,w,same_nr_types,N)

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

allscreen()

P1_nl=avg_samples(q_nl,create_timecell(ro,length(q_nl)));
P2_nl=avg_samples(p_nl,create_timecell(ro,length(p_nl)));


P1=avg_samples(q,create_timecell(ro,length(q)));
P2=avg_samples(p,create_timecell(ro,length(p)));

subplot(6,2,1)
plot(cell2mat(create_timecell(ro,1)),P2_nl(1,:))
title('Baseline HPC')
win1=[min(P2_nl(1,:)) max(P2_nl(1,:)) min(P2(1,:)) max(P2(1,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
subplot(6,2,3)
plot(cell2mat(create_timecell(ro,1)),P1_nl(1,:))
win1=[min(P1_nl(1,:)) max(P1_nl(1,:)) min(P1(1,:)) max(P1(1,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])

subplot(6,2,5)
plot(cell2mat(create_timecell(ro,1)),P2_nl(2,:))
win1=[min(P2_nl(2,:)) max(P2_nl(2,:)) min(P2(2,:)) max(P2(2,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
title('Baseline PFC')
subplot(6,2,7)
plot(cell2mat(create_timecell(ro,1)),P1_nl(2,:))
win1=[min(P1_nl(2,:)) max(P1_nl(2,:)) min(P1(2,:)) max(P1(2,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])

subplot(6,2,9)
plot(cell2mat(create_timecell(ro,1)),P2_nl(3,:))
win1=[min(P2_nl(3,:)) max(P2_nl(3,:)) min(P2(3,:)) max(P2(3,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
title('Baseline PAR')
subplot(6,2,11)
plot(cell2mat(create_timecell(ro,1)),P1_nl(3,:))
win1=[min(P1_nl(3,:)) max(P1_nl(3,:)) min(P1(3,:)) max(P1(3,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])




subplot(6,2,2)
plot(cell2mat(create_timecell(ro,1)),P2(1,:))
title('Plusmaze HPC')
win1=[min(P2_nl(1,:)) max(P2_nl(1,:)) min(P2(1,:)) max(P2(1,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
subplot(6,2,4)
plot(cell2mat(create_timecell(ro,1)),P1(1,:))
win1=[min(P1_nl(1,:)) max(P1_nl(1,:)) min(P1(1,:)) max(P1(1,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])

subplot(6,2,6)
plot(cell2mat(create_timecell(ro,1)),P2(2,:))
win1=[min(P2_nl(2,:)) max(P2_nl(2,:)) min(P2(2,:)) max(P2(2,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
title('Plusmaze PFC')
subplot(6,2,8)
plot(cell2mat(create_timecell(ro,1)),P1(2,:))
win1=[min(P1_nl(2,:)) max(P1_nl(2,:)) min(P1(2,:)) max(P1(2,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])

subplot(6,2,10)
plot(cell2mat(create_timecell(ro,1)),P2(3,:))
title('Plusmaze PAR')
win1=[min(P2_nl(3,:)) max(P2_nl(3,:)) min(P2(3,:)) max(P2(3,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])
subplot(6,2,12)
plot(cell2mat(create_timecell(ro,1)),P1(3,:))
win1=[min(P1_nl(3,:)) max(P1_nl(3,:)) min(P1(3,:)) max(P1(3,:))];
win1=[(min(win1)) (max(win1))];
ylim(win1)
xlabel('Time (s)')
ylabel('uV')
xlim([-.1 .1])


end
