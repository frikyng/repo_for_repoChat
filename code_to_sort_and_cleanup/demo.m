obj.rendering = false;


figure(666);clf();set(gcf,'Color','w');
t = obj.timescale.global_timescale;
valid = find(~cellfun(@isempty, obj.spiketrains.spike_estimate),1,'first');

ax1 = subplot(5,1,1);
plot(t, obj.binned_data.median_traces);hold on;title('Median per bin per group, and inferred spikes for bin 1 (red)'); hold on;
bsl = nanmin(nanmin(obj.binned_data.median_traces));
R = nanmax(range(obj.binned_data.median_traces)/20);
plot_raster(obj.spiketrains.spike_estimate{valid}, ones(size(obj.spiketrains.spike_estimate{valid}))*bsl - R, R,'r',2);

ax2 = subplot(5,1,2);
R = max(range(obj.binned_data.median_traces));
norm_vmr = obj.variability.index_of_disp;
plot(obj.event_fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
hold on;plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');

ax3 = subplot(5,1,3);
weighted_averages = obj.get_weight_map();title('Cell weighted average per component'); hold on;
plot(obj.timescale.global_timescale, weighted_averages);legend();

ax4 = subplot(5,1,4);
hold on;plot(obj.timescale.global_timescale, obj.variability.corr_results,'Color',[0.9,0.9,0.9]);
hold on;plot(obj.timescale.global_timescale, nanmean(obj.variability.corr_results, 2),'r');title('Pairwise temporal correlation between traces');ylim([-1.1,1.1]);
plot([0,nanmax(obj.timescale.global_timescale)],[-1,-1],'k--');plot([0,nanmax(obj.timescale.global_timescale)],[1,1],'k--')

ax5 = subplot(5,1,5);
[bouts, behav] = obj.get_activity_bout('encoder'); starts = bouts(1:2:end); stops = bouts(2:2:end);
plot(t, behav);title('Running speed and active bouts & inferred spiketrain(red)'); hold on;
plot_raster(obj.spiketrains.spike_estimate{valid}, ones(size(obj.spiketrains.spike_estimate{valid}))- 5, 5,'r',2);
for el = 1:numel(starts)
    x = t([starts(el),starts(el),stops(el),stops(el)]);
    y = [0,nanmax(behav),nanmax(behav),0];
    patch('XData',x,'YData',y,'FaceColor','red','EdgeColor','none','FaceAlpha',.1);hold on
end
xlabel('Recording time (s)')

linkaxes([ax1, ax2, ax3, ax4, ax5],'x')


function h = plot_raster(x,y,width,color,linew)
    if nargin<1; help plot_raster; return; end
    if nargin<2 || isempty(y); y=x; x=1:length(y); end
    if nargin<3 || isempty(width); width=1; end
    if nargin<4 || isempty(color); color='k'; end
    if nargin<5 || isempty(linew); linew=1; end
    if length(x)==1 && length(y)>1
        x=ones(size(y))*x;
    end
    if length(y)==1 && length(x)>1
        y=ones(size(x))*y;
    end
    plot([x(:)'; x(:)'], [y(:)'+width/2; y(:)'-width/2],'-','color',color,'linewidth',linew);
end