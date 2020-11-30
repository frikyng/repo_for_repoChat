function assess_variability(all_pks, all_concat, peak_times, global_timescale)
    sr = nanmedian(diff(global_timescale));
    vmr = nanvar(all_pks,[],2)./nanmean(all_pks, 2);
    cv  = nanstd(all_pks,[],2)./nanmean(all_pks, 2); % (maybe chack snr at one point?  mu / sigma)
    fano = []; % windowed VMR. usually for spike trains
    [~, idx] = sort(vmr,'descend');
    figure(123); cla();ylim([0, nanmax(all_pks(:))]); hold on;
    
    %     for event = idx'
    %         show_event(all_concat, round(peak_times/sr), event);
    %         drawnow;%pause(0.1)
    %     end
    
    figure(1009);cla();plot(all_pks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')
    
    R = max(range(all_concat));
    %norm_vmr = vmr/range(vmr);
    norm_vmr = vmr/mean(vmr);
    figure(1010);clf();
    ax1 = subplot(2,1,1);plot(global_timescale, all_concat); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
    ax2 = subplot(2,1,2);plot(peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
    plot([global_timescale(1), global_timescale(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
    plot([global_timescale(1), global_timescale(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
    hold on;plot([global_timescale(1), global_timescale(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');
    linkaxes([ax1, ax2], 'x');
    
    figure();histogram(norm_vmr, round(10*range(norm_vmr)/std(norm_vmr)))
    
    figure(1011);cla();scatter(nanmedian(all_pks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
end

