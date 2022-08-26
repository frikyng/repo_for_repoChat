
function rise = get_rise_time(exp)
test = diff(smoothdata(exp.binned_data.median_traces,'gaussian',[10,0]));
bsl = rms(test,'OmitNan')*2;
test(abs(test) < bsl) = NaN;
rise = test > 0;
decay = test < 0;
rise = double(any(rise, 2));
%rise(rise == 1) = -10;
%rise(rise == 0) = NaN;
% figure();plot(exp.binned_data.global_median);hold on;scatter(1:numel(rise),rise )
end