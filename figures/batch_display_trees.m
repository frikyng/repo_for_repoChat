%% Displays all tree next to each other using analyzed data

function out = batch_display_trees(top_folder, settings_file_path, meta_analysis_path, mode)
    if nargin < 2 || isempty(settings_file_path) 
        [~,~,~,~,top_folder] = parse_paths(top_folder, true);
        settings_file_path = [top_folder, 'settings.txt'];
    end
    if nargin < 4 || isempty(mode) 
        mode = 'corr';
    end
    load(meta_analysis_path);%"C:\Users\vanto\Desktop\new_meta_june_2020_v6\summary.mat")
    
    %filter = 'whatever'
    %plot_group_correlation_data(results, filter);
    
    %% Load batch processing settings
    batch_params = get_paths_from_settings(settings_file_path); 

    %% Now generate the figures
    all_trees   = {};
    all_soma    = {};
    all_values  = {};
    
    extracted = find(~cellfun(@isempty, results.dimensionality));
    for expe = extracted
        if strcmp(mode, 'corr')
            [all_trees{expe}, all_soma{expe}, all_values{expe}] = results.plot_corr_tree(expe);
        elseif strcmp(mode, 'dim')
            [all_trees{expe}, all_soma{expe}, all_values{expe}] = results.plot_dim_tree(0, expe); 
        elseif strcmp(mode, 'dist')
            [all_trees{expe}, all_soma{expe}, all_values{expe}] = results.plot_dist_tree(0, expe); 
        elseif strcmp(mode, 'order')
            [all_trees{expe}, all_soma{expe}, all_values{expe}] = results.plot_ord_tree(0, expe); 
        end
        close all;
    end
    
    figure(); axis equal;set(gcf,'Color','w');
    offset = 0;
    col = lines(numel(all_values));
    row = 0;
    max_w = 5000;
    for expe = 1:numel(all_values)
        if ~isempty(all_values{expe})
            if offset > (max_w-500) || expe == 1
                offset = 0;
                if expe > 1
                    row = row + 700;
                end
    %             plot([0,max_w],[row,row],'k-'); hold on;
    %             plot([0,max_w],[row+180,row+180],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
    %             plot([0,max_w],[row+360,row+360],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
    %             plot([0,max_w],[row+600,row+600],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
            end
            shifted_x_coor = all_values{expe}{1}(1,:); 
            shifted_z_coor = all_values{expe}{3}(1,:);
            zero_shifted_coor = -min(shifted_x_coor) + offset;
            shifted_x_coor = shifted_x_coor + zero_shifted_coor;
            X = shifted_x_coor;
            Y = shifted_z_coor+row;
            %Z = zeros(size(X));
            Z = all_values{expe}{2}(1,:);
            S = all_values{expe}{4}(1,:);
            
            hold on;
            HP = surface('XData'    ,[X;X],...
                         'YData'    ,[Y;Y],...
                         'ZData'    ,[Z;Z],...
                         'CData'    ,[S;S],...
                         'facecol'  ,'no',...
                         'edgecol'  ,'interp',...
                         'linew'    ,1); hold on
                     
            scatter(all_soma{expe}(1) + zero_shifted_coor,all_soma{expe}(3)+row,50,'k','filled'); hold on;
            offset = max(shifted_x_coor);
        end
    end
    ylim([0, row+800]);
    xlim([0,max_w]); hold on;set(gca, 'Ydir', 'reverse');colorbar                                                       
end

function plot_group_correlation_data(results, filter)
            % Plot Group pearson  
            
%     if ~isempty(filter)
%         filter = cellfun(@(x) size(x, 2), [results.cumsum]) < 9;
%         results = structfun(@(x) x(filter), results, 'UniformOutput', false);
%     end       
%             
            
    max_gp = max(cellfun(@(x) size(x,2) ,results.crosscorr));
    results.crosscorr = results.crosscorr(~cellfun(@isempty, results.crosscorr));
    results.cumsum = results.cumsum(~cellfun(@isempty, results.cumsum));
    results.peaks = results.peaks(~cellfun(@isempty, results.peaks));
    

    
    
    all_cc = cellfun(@(x) [x(1:size(x,1),1:size(x,2)), NaN(size(x,1),max_gp-size(x,2))], results.crosscorr, 'UniformOutput', false);
    all_cc = cellfun(@(x) [x(1:size(x,1),1:size(x,2)); NaN(max_gp-size(x,1),size(x,2))], all_cc, 'UniformOutput', false);
    all_cc = cat(3,all_cc{:});    
    figure();imagesc(nanmean(all_cc,3));axis equal;colorbar
    
    figure();imagesc(sum(~isnan(all_cc),3));axis equal;colorbar
    
    
    first_gp_cc = cellfun(@(x) x(2,:), results.crosscorr, 'UniformOutput', false);
    first_gp_cc = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,max_gp-size(x,2))], first_gp_cc, 'UniformOutput', false);
    first_gp_cc = vertcat(first_gp_cc{:});
    figure();bar(nanmean(first_gp_cc));hold on; errorbar(nanmean(first_gp_cc), nanstd(first_gp_cc)./sqrt(sum(~isnan(first_gp_cc))), 'k') ;

    % Std per event (need behaviour) 
    %test = cellfun(@(x) nanstd(x'), results.peaks, 'UniformOutput', false); % std per event

    % Variability per group (from mean)
    test = cellfun(@(x) nanvar(x)./nanmean(x), results.peaks, 'UniformOutput', false); % variance from mean, per subgroup
    test = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,max_gp-size(x,2))], test, 'UniformOutput', false);
    test = vertcat(test{:});
    figure();bar(nanmean(test));hold on; errorbar(nanmean(test), nanstd(test)./sqrt(sum(~isnan(test))), 'k') 
end