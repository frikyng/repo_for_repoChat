%% Apply wavelet denoising.

function [data, missing] = wavelet_denoise(data)
    missing             = isnan(data);
    data                = fillmissing(data,'linear');
    skip                = ~all(isnan(data)); % When one trace is full of NaN, we can't analyse it
    data(:,skip)        = single(wdenoise(double(data(:,skip))));
    data(missing)       = NaN;
end

