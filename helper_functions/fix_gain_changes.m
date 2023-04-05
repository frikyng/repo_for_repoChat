
%% Rescale traces
% - if expe.detrend == -1, for each ROI, we rescale the F0 of each trial to the median F0 of the ROI
% - if expe.detrend > 0  , for each ROI, we compute a fit of the baseline
%       drift using a polynomial fit of order "expe.detrend£ (e.g.
%       expe.detrend = 1 is a linear fit). This is then converted into a
%       "gain function" that will amplify or dim the signal accordingly
% - if expe.detrend == 'auto', breakpoints are identified automatically,
%   then a linear detrending per block is applied
% - if breakpoints_idx is not empty, the fit is done per block of trials,
%   as defined by the breakpoints_idx limits.

function extracted_traces = fix_gain_changes(expe, extracted_traces)
    temp        = vertcat(extracted_traces{:}); % get traces
    real_low    = nanmin(temp(:)); % real dark noise level (won't change if gain change as it is not amplified by PMT's)
    temp        = temp - real_low + eps; % remove dark noise. eps to prevent division by 0
    max_v       = prctile(temp(:),90);%max(cellfun(@(x) nanmax(x(:)), extracted_traces));
    
    detrending_plot = expe.rendering;
    detrending_plot = false;
    if detrending_plot
        figure(667788);clf();ax1 = subplot(1,3,1);imagesc(temp'); caxis([0, max_v]);
    end
    
    if ischar(expe.detrend) && strcmp(expe.detrend, 'auto')
        F0         = cell2mat(cellfun(@(x) prctile(x - real_low + eps, 1)', extracted_traces, 'UniformOutput', false));    
        ref        = nanmedian(F0, 2);
        scaling_fact = F0 .\ ref;    
        [expe.breakpoints,~] = findchangepts(scaling_fact(~all(isnan(scaling_fact),2),:));
        expe.detrend = 1;
    end
    
%     F0         = cell2mat(cellfun(@(x) prctile(x - real_low + eps, 1)', extracted_traces, 'UniformOutput', false));    
%     ref        = nanmedian(F0, 2);
%     scaling_fact = F0 .\ ref;    
%     invalid = range(scaling_fact,2) > 0.2;
%     expe.bad_ROI_list = find(invalid);
%     
%     
    if expe.detrend == -1   
        %% Detrending using linear regression
        F0         = cell2mat(cellfun(@(x) prctile(x - real_low + eps, 1)', extracted_traces, 'UniformOutput', false));    
        ref        = nanmedian(F0, 2);
        scaling_fact = F0 .\ ref;    
        for trial = 1:numel(extracted_traces)
            for row = 1:size(extracted_traces{1},2)                
                extracted_traces{trial}(:, row) = (extracted_traces{trial}(:, row) - real_low + eps).*scaling_fact(row,trial) + real_low - eps;
            end
            slope_gain{1}{trial} = scaling_fact(:, trial)';
        end
    elseif expe.detrend == -2   
        for trial = 1:numel(extracted_traces)
            [extracted_traces{trial}, slope_gain{1}{trial}] = normalize(extracted_traces{trial}, 1, 'medianiqr');
        end
    elseif expe.detrend == -3 

        %% dF/F0 for the signal channel
        Fb                          = repmat(nanmin(temp(:)),size(temp,1),1)-1;  
        temp                        = temp - Fb;
        F0                          = movmin(temp, 300,1) + eps;%repmat(signal_baseline,1,timepoints);
        temp                        = (temp - F0) ./ (F0);
        
        counter = 1;
        for trial = 1:numel(extracted_traces)
            tp = expe.timescale.tp(trial);
            extracted_traces{trial} = temp(counter:counter+tp-1,:);
            counter = counter + tp;
        end
        return
        
    else
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
    
    if detrending_plot
        ax2 = subplot(1,3,2);imagesc(cat(1, extracted_traces{:})');%hold on;caxis([-max_v, max_v]);
        ax3 = subplot(1,3,3);imagesc(temp' - cat(1, extracted_traces{:})');hold on;caxis([0, max_v]); 
        linkaxes([ax1, ax2, ax3], 'xy')  
        if expe.detrend
            figure(667799);clf();
            if expe.detrend == -1  
                subplot(2,1,1);imagesc(vertcat(slope_gain{1}{:})');ylabel('ROI');xlabel('trial');title('Detrending - linear gain per trial');
            elseif expe.detrend == -2  
                subplot(2,1,1);imagesc(vertcat(slope_gain{1}{:})');ylabel('ROI');xlabel('trial');title('Detrending - median interquartile normalization');        
            else
                all_gains = cellfun(@(x) horzcat(x{:}), slope_gain, 'UniformOutput', false);
                subplot(2,1,1);plot(horzcat(all_gains{:}));ylabel('gain');xlabel('timpoints');title('Detrending - gain correction');
            end
            hold on; subplot(2,1,2);plot(nanmedian(temp,2),'r'); hold on; plot(nanmedian(cat(1, extracted_traces{:}),2), 'k')
        end
    end
end

