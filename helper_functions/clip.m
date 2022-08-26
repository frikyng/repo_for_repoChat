function in = clip(in,thr)
    if nargin < 2 || isempty(thr)
        thr = [0, 1];
    end
    
    if numel(thr) == 1
        in(in > thr(1)) = thr(1);
    else
        in(in < thr(1)) = thr(1);
        in(in > thr(2)) = thr(2);
    end
end