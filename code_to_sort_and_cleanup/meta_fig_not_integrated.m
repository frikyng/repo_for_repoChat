
for expe_idx = 25%53%78%80:94
    data = a.experiments(expe_idx).rescaled_traces;
    
    %data = vertcat(a.experiments(expe_idx).rescaled_traces);
    
    t = a.experiments(expe_idx).event_fitting.peak_pos;
    mean_trace = nanmedian(data, 2);
    %t = t(mean_trace(t) > 50);
    amp = NaN(size(t,2), size(data,2));
    chance_amp = NaN(size(t,2), size(data,2));
    data_sm = smoothdata(data, 'gaussian', 10);
    for trace_idx = 1:size(data, 2)
        amp(:, trace_idx) = data_sm(t, trace_idx); 
        chance_amp(:, trace_idx) = data_sm(randi(size(data, 1), [1, numel(t)]), trace_idx);  
    end


    usable = ~all(isnan(data),1);
    test_sample = data;
    noise_level = NaN(1,size(data, 2));
    %test_sample = data(~isnan(mean(data(:,usable),2)),usable); %data no nan
    for trace_idx = find(usable)
        sm = smoothdata(test_sample(:, trace_idx), 'gaussian', 10);
        sm(sm > rms(sm) | isnan(sm)) = [];
        noise_level(trace_idx) = rms(sm)*2;
    end
    %test_sample = test_sample(~isnan(mean(test_sample,2)),:); %data no nan
    %noise_level = rms(test_sample);

    mean_amp = nanmedian(data(t,:), 2);
    for trace_idx = 1:size(data, 2)
        [pk_amp, pk_loc] = findpeaks(smoothdata(double(data(:,trace_idx)), 'gaussian', [20,0]), 'MinPeakProminence' ,noise_level(trace_idx));
%         figure(123);cla();findpeaks(smoothdata(double(data(:,trace_idx)), 'gaussian', [20,0]), 'MinPeakProminence' ,noise_level(trace_idx));
%         pause(0.5)
    end
    
    figure(125);clf();
    subplot(2,2,1);bar(sum(amp > noise_level*2) ./ numel(t) *100); ylabel('% of sucCess');xlabel('ROI #');title('% of events above noise')
    subplot(2,2,2);bar(sum(chance_amp > noise_level*2) ./ numel(t) *100); ylabel('% of sucCess');xlabel('ROI #');title('Chance level')
    subplot(2,2,[3,4]);plot(mean_trace)
    

    figure(126);clf();
    subplot(2,3,1);c1 = corrcoef(amp);imagesc(c1);caxis([0,1]);axis image;title('Correlation of peaks amp');
    subplot(2,3,2);c2 = corrcoef(data,'Row','pairwise');imagesc(c2);caxis([0,1]);axis image;title('Correlation of traces')
    subplot(2,3,3);imagesc(c1-c2);axis image;caxis([0,1]);colormap(redblue)
    s2 = subplot(2,3,4);a.experiments(expe_idx).ref.plot_value_tree(c1(:,1),'','','','',s2);caxis([0,1]);
    s3 = subplot(2,3,5);a.experiments(expe_idx).ref.plot_value_tree(c2(:,1),'','','','',s3);caxis([0,1]);
    s = subplot(2,3,6);a.experiments(expe_idx).ref.plot_value_tree(nanmean(amp),'','','','',s);colormap(viridis)
 

    
    
    pause(1)
end