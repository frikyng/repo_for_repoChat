%% Apply wavelet denoising.

function [data, missing] = wavelet_denoise(data, varargin)
    missing             = isnan(data);
    data                = fillmissing(data,'linear');
    keep                = ~all(isnan(data)); % When one trace is full of NaN, we can't analyse it
    data(:,keep)        = single(wdenoise(double(data(:,keep)), varargin{:}));
    data(missing)       = NaN;
end

%     figure();hold on;
%     plot(data,'k');hold on;
%     plot(wdenoise(double(data)),'r');hold on
%     plot(wdenoise(double(data),'DenoisingMethod','BlockJS'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','FDR'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','Minimax'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','Sure'));hold on
%     legend({'original','Bayes','BlockJS','FDR','Minimax','Sure'}) ; % BAYES, FDR and SURE look good