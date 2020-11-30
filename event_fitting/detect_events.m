function [pk_val, pk_loc, pk_t, widths_in_pt, denoised] = detect_events(data, timescale, amp_factor)
    if nargin < 2 || isempty(timescale)
        timescale = 1:size(data, 1);
        tlabel = '(frames)';
    else
        tlabel = '(s)';
    end
    if nargin < 3 || isempty(amp_factor)
        amp_factor = 1;
    end
    
    %% Filter signal
    missing             = isnan(data);
    data                = fillmissing(data,'linear');
    denoised            = wdenoise(double(data));
    bsl                 = data - denoised;
    %thr                 = rms(bsl) * thr_factor;
    %figure(25);cla();envelope(bsl,20,'peak');
    [a, b] = envelope(bsl,20,'peak');
    thr = (prctile(a-b,99) - prctile(a-b,1))*3*amp_factor;
    fprintf('PLEASE REVIEW THRESHOLDING CRITERIA BEFORE FINAL ANALYSIS\n')
    
    %figure();hold on;
    %plot(data,'k');hold on;
    %plot(wdenoise(double(data)),'r');hold on
    %plot(wdenoise(double(data),'DenoisingMethod','BlockJS'));hold on
    %plot(wdenoise(double(data),'DenoisingMethod','FDR'));hold on
    %plot(wdenoise(double(data),'DenoisingMethod','Minimax'));hold on
    %plot(wdenoise(double(data),'DenoisingMethod','Sure'));hold on
    %legend({'original','Bayes','BlockJS','FDR','Minimax','Sure'}) ; % BAYES, FDR and SURE look good
    
    %% Find peaks per group (if multiple ones)
    denoised(missing)   = NaN;
    pk_val              = {};
    pk_loc              = {};
    pk_t                = {};
    widths_in_pt              = cell(size(denoised, 2), 1);

    for gp = 1:size(denoised, 2)
        figure(1002);cla();hold on;
        [pk_val{gp}, pk_loc{gp}, widths_in_pt{gp}] = findpeaks(denoised(:, gp),'MinPeakProminence',thr(gp),'WidthReference','halfheight'); % 100 ms min width
        pk_t{gp} = timescale(pk_loc{gp});
        plot(timescale, data(:, gp), 'k'); hold on;
        plot(timescale, denoised(:, gp), 'r'); hold on;
        scatter(pk_t{gp}, pk_val{gp}, 'kv', 'filled');
        ylabel('Amplitude')
        xlabel(['Time ',tlabel]);
        title('Event detection')
    end
end

