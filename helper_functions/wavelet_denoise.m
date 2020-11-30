function data = wavelet_denoise(data)
    missing             = isnan(data);
    data                = fillmissing(data,'linear');
    skip                = ~all(isnan(data)); % QQ not sure whene we have full NaN series, but we can't denoise them with wdenoise
    data(:,skip)        = single(wdenoise(double(data(:,skip))));
    data(missing)       = NaN;
end

