% demo simple : 27 ; 2019-09-24_exp_1
% demo regulier: 37 ; 2019-10-03_exp_5  (or 9, 2019-09-05_exp_1
% demo super fast events: 15 ; 2019-09-10_exp_1
% demo fast and slow kinetics : 11 : 2019-09-07_exp_7 (alos interresting 2019-10-08_exp_2)
% demo burst : 13, 2019-09-09_exp_2
% horrible with bursts and low sampling rate: 5 ; 2019-08-08_exp_1 
% TOFIX folder 3, 2019-08-02_exp_1 had an initial extraction issue, where apical dendrite data were not even saved. go back to the extraction process
% SCALING ISSUE 2019-10-04_exp_1
% and maybe 2019-10-07_exp_1 too


%% Set some analysis options
opengl software
save_data   = true;
demo        = 0; %0 / false for no demo ; 1 for fitting demo ; 2 for fitting and scaling demo
filter      = '';
condition   = {'distance', 100};
    
%% Define the folders to use
%% source_folder points to the output of batch_process_ribbon_scan()
% setting_file_path = 'D:\Curated Data\settings.txt';
% source_folder = 'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER\Ribbon Scan paper\code for paper\test7\meta\';
% export_folder = 'C:\Users\vanto\Desktop\new_meta_june_2020_v12';



if isvarname('results') && isa(results, 'arboral_scan_meta_analysis')
    
elseif isfile([export_folder, 'summary.mat'])
    %% Resume analysis
    load([export_folder, 'summary.mat']);
else
    %% Start new analysis
    results                 = arboral_scan_meta_analysis(source_folder, export_folder);
    setting_file_path = SETTINGS_FILE%'D:\Curated Data\settings.txt';
    source_folder = [EXPORT_PATH, '\meta']%'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER\Ribbon Scan paper\code for paper\26-10-2020\meta\';
    export_folder = [pwd, '\new_meta_nov_2020'];
    export_folder = parse_paths(export_folder);
end
data_folders_per_exp    = results.filter_expe_subset();



data_folders_per_exp = data_folders_per_exp(results.need_update)

failed_analysis         = {};
failed_factoran         = {};
whitebg('w')

%% Now, do the analysis expe-by-expe
for expe = 1:numel(data_folders_per_exp) % review expe 23 for fitting warning, saturated peak scaling and one group having an overestimated baseline
    original_expe_idx = expe; % for the saving part 
    
    %% Close existing figures
    close all
    fprintf(['Processing experiment #',num2str(expe),'\n'])
    
    %try
    %% Prepare fields some useful information
    expe = results.add_experiment(data_folders_per_exp{original_expe_idx});
        
    results.filter_win = [0, 100];
    
    %% Load and concatenate traces for the selected experiment
    results.load_extracted_data(expe);   % this also sets the current expe #
    
    
    %% Check for consistent number of ROIs
    quality_control(1, results.extracted_traces{expe})
    %quality_control(2, cell2mat(all_sr))

    %% Prepare binning of ROIs based on specific grouping condition
    [bins, metrics, bin_legend] = results.prepare_binning(condition);

    %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
    results.rescale_traces();
    figure();plot(nanmedian(results.get_rescaled_traces(),2))
    arrangefigures([1,2]); 

    %% Create median trace per bins
    results.set_median_traces()

    %% Plot the mean trace for each bin  
    figure(1001);cla();plot(results.timescale{expe}.global_timescale, results.binned_data{expe}.median_traces); hold on;legend(results.binned_data{expe}.bin_legend); hold on;title('mean scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
    arrangefigures([1,2]); 
    
    %% Get summary covariance plots
    results.similarity_plot();
    arrangefigures([1,2]); 

    %% Correct for decay to avoid overstimating peak amplitude
    if isempty(results.event_fitting{expe})
        results.event_fitting{expe} = detect_and_fit_events(results.binned_data{expe}.median_traces, results.timescale{expe}.global_timescale, demo, results.binned_data{expe}.bin_legend);
        arrangefigures([1,2]);
    end

    % peak_times = results.event_fitting{expe}.peak_times
    % peak_pos = results.event_fitting{expe}.peak_pos
    % corrected_amp = results.event_fitting{expe}.post_correction_peaks
    
    %% Detect and display peak histogram distribution (mean and individual groups)
    results.plot_events_distribution();
    
    %% Study bAP heterogeneity
    results.assess_variability()

    %% Check how correlated are the peaks between different parts of the tree
    if isempty(results.crosscorr{expe})
        results.crosscorr{expe} = corrcoef([nanmean(results.event_fitting{expe}.post_correction_peaks, 2), results.event_fitting{expe}.post_correction_peaks]);
    end
    results.plot_cc(); arrangefigures([1,2]); 

    %% Plot a tree correlation metric     
    %ROIs_list  = unique([all_ROI_ids_per_bin{:}]);
    %[fixed_tree, soma_location, ROIs_list, listing, batch_params] = rebuild_tree(data_folders_per_exp{original_expe_idx}(1), setting_file_path, true);
    
    %% Assign value of group to these ROIs
    results.plot_corr_tree(); arrangefigures([1,2]);

    cross_validate = false;
    results.get_dimensionality(cross_validate);
    
    %% Plot weight-tree for each component
    for comp = 1:5
        results.plot_dim_tree(comp);
    end
    
    %% Plot a map of tree weights by ROI number for each component
    results.plot_weight_map();

    %% Plot strongest componenet per ROI
    results.plot_strongest_comp_tree();

    arrangefigures([1,2]);
    
    %% Optionally, if external variables need an update
    %results.update_external_metrics(60)

    %% Save figure
    if save_data
        p = get(groot,'DefaultFigurePosition');
        folder = [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/'];
        if isfolder(folder)
            rmdir(folder,'s');
        end
        mkdir(folder);
        for fig_idx = [1001:1021, 1081, 1082, (10020+1):(10020+size(results.event_fitting{expe}.events, 1)), (20020+1):(20020+5)]
            f = figure(fig_idx);
            set(f, 'Position', p)
            try
                saveas(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                saveas(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                savefig(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
            catch % for figures with subplots
                try
                    saveas(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                    saveas(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                    savefig(f, [results.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
                end
            end
        end
        name = parse_paths([folder,'/',f.Tag,data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.mat']);
        batch_params = results.general_info{expe}.arboreal_scan.batch_params;
        
        %% Save every time in case of failure.
        save(parse_paths([results.export_folder, 'summary.mat']), 'results','-v7.3')
        %save(name, 'results.binned_data{expe}.median_traces','all_fits', 'all_pks', 'peak_times','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','-v7.3');   
        %save(name, 'all_traces_per_bin','all_traces_per_bin','all_traces_concat','all_fits', 'all_pks', 'peak_times','cc','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','ROIs_list','dimensionality','-v7.3');   
    end

    %catch
    %   failed{expe} = data_folders_per_exp{original_expe_idx}; 
	%end
    close all
end

if save_data
    save(parse_paths([results.export_folder, 'summary.mat']), 'results','-v7.3')
end

% Plot Group pearson  
max_gp = max(cellfun(@(x) size(x,2) ,results.crosscorr));
results.crosscorr = results.crosscorr(~cellfun(@isempty, results.crosscorr));
% results.cumsum = results.cumsum(~cellfun(@isempty, results.cumsum));
% results.peaks = results.peaks(~cellfun(@isempty, results.peaks));
test = cellfun(@(x) x(1,:), results.crosscorr, 'UniformOutput', false);
test = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,max_gp-size(x,2))], test, 'UniformOutput', false);
test = vertcat(test{:});
figure();bar(nanmean(test));hold on; errorbar(nanmean(test), nanstd(test), 'k') ;

% Std per event (need behaviour) 
test = cellfun(@(x) nanstd(x'), results.peaks, 'UniformOutput', false); % std per event

% Variability per group (from mean)
test = cellfun(@(x) nanstd(x.post_correction_peaks - nanmean(x.post_correction_peaks, 2)), results.event_fitting, 'UniformOutput', false); % variance from mean, per subgroup
test = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,max_gp-size(x,2))], test, 'UniformOutput', false);
test = vertcat(test{:});
figure();bar(nanmean(test));hold on; errorbar(nanmean(test), nanstd(test), 'k') 



