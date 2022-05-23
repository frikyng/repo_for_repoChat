function extracted_traces = fix_gain_changes(expe, extracted_traces)
    temp        = vertcat(extracted_traces{:}); % get traces
    real_low    = 0;%nanmin(temp(:)); % real dark noise level (won't change if gain change)
    temp        = temp - real_low; % remove dark noise
    tp_start    = cumsum([1, expe.timescale.tp]); % pt start of trials
    breakpoints_idx = find(cellfun(@(x) ~isempty(expe.breakpoints) && contains(x, expe.breakpoints),expe.updated_data_path));
    blocks      = [1, breakpoints_idx, numel(extracted_traces)+1]; % pt range of the blocks
    pts         = [1,tp_start(breakpoints_idx),size(temp,1)+1]; % pt start of blocks
    bsl_value   = [];
    slope_gain  = {};

    for block_idx = 1:(numel(pts)-1)
        tmp_traces = temp(pts(block_idx):pts(block_idx+1)-1,:);

        %% Fix slop first since we'll fix the gain later
        local_tp        = size(tmp_traces,1);

        if expe.detrend
            if any(isnan(tmp_traces(:)))
                tmp_traces  = fillmissing(tmp_traces,'movmean',5);
                to_use      = ~all(isnan(tmp_traces), 1);
            else
                to_use      = true(1,size(tmp_traces, 2));
            end

            %                     %% Get slope for current trace
            %                     p = polyfit(repmat(1:local_tp,sum(to_use),1),movmin(tmp_traces(:,to_use),50)',1);
            %                     slope = (local_tp*p(1)) ./ nanmean(tmp_traces(:)); % slope in percent
            %                     slope_gain{block_idx} = linspace(-slope/2,slope/2,local_tp); 

            %                     x = 1:local_tp;
            %                     a = nanmean(movmin(tmp_traces(:,to_use),50),2);
            %                     out = fit(x',a,'a + b*log(c - x)','Lower',[0,1,numel(a)],'Upper',[max(a),100,numel(a)*2]);toc
            %                     slope_gain{block_idx} = out.a + out.b*log(out.c-x);
            %                     %slope_gain{block_idx}(slope_gain{block_idx} < nanmin(a)) = nanmin(a);
            %                     slope_gain{block_idx} = slope_gain{block_idx} ./ mean(slope_gain{block_idx}) - 1;;

            if islogical(expe.detrend)
                poly_order = 1;
            else
                poly_order = expe.detrend;
            end
            x = 1:local_tp;
            p = polyfit(repmat(x,sum(to_use),1),movmin(tmp_traces(:,to_use),[300, 0])',poly_order);
            % figure(555);cla();plot(nanmean(movmin(tmp_traces(:,to_use),300)')); hold on;plot(x, polyval(p,x))
            slope_gain{block_idx} = polyval(p,x);
            slope_gain{block_idx} = slope_gain{block_idx} ./ nanmean(slope_gain{block_idx}) - 1;

            %                     out = fit(x',a,'a + b*log(c - x)','Lower',[0,1,numel(a)],'Upper',[max(a),100,numel(a)*2]);toc
            %                     slope_gain{block_idx} = out.a + out.b*log(out.c-x);
            %                     %slope_gain{block_idx}(slope_gain{block_idx} < nanmin(a)) = nanmin(a);
            %                     slope_gain{block_idx} = slope_gain{block_idx} ./ mean(slope_gain{block_idx}) - 1;;

            %% Temporarily correct the block
            tmp_traces = tmp_traces - tmp_traces.* slope_gain{block_idx}';

            %% Reslice per trial
            slope_gain{block_idx} = mat2cell(slope_gain{block_idx},1, expe.timescale.tp(blocks(block_idx):(blocks(block_idx+1)-1)));
        end
        if ~isempty(breakpoints_idx)
            bsl_value = [bsl_value; prctile(tmp_traces,1)];
        end
    end

    %% Get scaling if gain fix is required
    if ~isempty(breakpoints_idx)
        scaling     = bsl_value .\ bsl_value(1,:); % normalize to first group;
    end

    %% Now apply correction per trial
    for block_idx = 1:(numel(pts)-1)
        trial_idx                     = (blocks(block_idx)):(blocks(block_idx+1)-1);
        if expe.detrend
            extracted_traces(trial_idx)   = cellfun(@(x, y) x - x.*y', extracted_traces(trial_idx),slope_gain{block_idx}, 'uni',false);
        end
        if ~isempty(breakpoints_idx)
            extracted_traces(trial_idx)   = cellfun(@(x) (x.*scaling(block_idx,:)), extracted_traces(trial_idx), 'uni',false);
        end
    end
end

