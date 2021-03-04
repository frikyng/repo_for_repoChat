%% Once you extracted data per experiment in 
opengl software
save_data = true;
demo        = 1; %0 / false for no demo ; 1 for fitting demo ; 2 for fitting and scaling demo

[export_folder, data_folders_per_exp, setting_file_path] = prepare_meta_analysis('C:\Users\vanto\Desktop\new_meta_june_2020_v5');

%% Now, do the analysis expe-by-expe
results = {};
failed = {};
whitebg('w')
for expe = 1:numel(data_folders_per_exp)

    %% Reload previous analysis (obtained with ribbon_spatial_analysis())
    [all_traces_per_bin, all_ROI_ids_per_bin, blacklist] = load_analyses(data_folders_per_exp{expe});
    
    %% TO INVESTIGATE : some recordings have less groups
    difference_of_bin  = nanmax(cellfun(@(x) numel(x), all_traces_per_bin)) - cellfun(@(x) numel(x), all_traces_per_bin); % get number of missing bins (because of different masks i guess)
    if any(difference_of_bin)
        error_box('Some odd analysis here... you should check that')
        for data_f_idx = 1:numel(all_traces_per_bin)
            if difference_of_bin(data_f_idx)
                complement                      = cell(1, difference_of_bin(data_f_idx));
                complement                      = cellfun(@(x) single(NaN(size(all_traces_per_bin{data_f_idx}{1},1), 1)), complement, 'UniformOutput', false);
                all_traces_per_bin{data_f_idx}  = [all_traces_per_bin{data_f_idx}, complement];
            end
        end
    end
    
    pb = cellfun(@isempty, all_traces_per_bin{1});
    if any(pb)
        fprintf([data_folders_per_exp{expe}(1).name,'\n'])
        all_traces_per_bin{1}
        for rec = 1:numel(all_traces_per_bin)
            ref = size(all_traces_per_bin{rec}{find(~pb, 1, 'first')});
            all_traces_per_bin{rec}(pb) = repmat({single(NaN(ref(1),1))}, 1, sum(pb));
        end
    end
    
    %% Reload data for the biggest series
    [~, maxloc]         = nanmin(difference_of_bin);
    source              = dir([data_folders_per_exp{expe}(maxloc).folder, '/', data_folders_per_exp{expe}(maxloc).name,'/*.mat']);
    data                = load([source.folder, '\', source.name]);
    sr                  = diff(data.params.timescale(1:2)); %% QQ INCORRECT
    tp                  = cellfun(@(x) size(x{1}, 1), all_traces_per_bin,'UniformOutput', false);  %% QQ INCORRECT
    trial_timescale     = cellfun(@(x) linspace(0, sr*x, x), tp', 'UniformOutput', false);
    global_timescale    = linspace(0, sr*sum(cell2mat(tp)), sum(cell2mat(tp)));
    
    %% Prepare some labels
    bin_legend          = cellstr(num2str([data.distance_range' , data.distance_range' + (data.distance_range(2) - data.distance_range(1))]));
    bin_legend          = strrep(bin_legend,'0  ','0 - ');

    %% Now find a scaling between the different groups
    [all_traces_per_bin, all_traces]  = scale_across_recordings(all_traces_per_bin, trial_timescale, bin_legend, demo);

    %% Regroup recordings   
    all_concat = vertcat(all_traces_per_bin{:});
    
    %% Adjust trace offset to have it well at 0;
    valid_gp            = find(~all(isnan(all_concat))); % NaN series happens when all tracs of a group are NaN, or if there are no traces for the group, but we can denoise them with wdenoise
    bsl                 = NaN(1, size(all_concat, 2));
    bsl(valid_gp)       = prctile(movmedian(all_concat(:, valid_gp), 20), 1) ;%prctile(wdenoise(double(fillmissing(all_concat(:, valid_gp),'linear')),'NoiseEstimate','LevelDependent'), 10); % uses aggressive denoising to get a proper baseline estimate
    all_concat          = all_concat - bsl;
    
    %% Plot the mean trace for each bin  
    figure(10010);cla();plot(global_timescale, all_concat); hold on;legend(bin_legend); hold on;title('mean scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
    
    %% Get summary covariance plots
    compute_similarity
    
    detect_events(nanmedian(all_concat, 2), global_timescale);

    %arrangefigures(1); 
    {expe, data_folders_per_exp{expe}(1).name}
    %pause(1)
end