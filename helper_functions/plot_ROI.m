function out = plot_ROI(roi)
    global all_traces data_folders mask_name
    hFig = figure(123);clf();plot(all_traces(:,roi));  hold on;  
    title(['ROI: ',num2str(roi)]);
    f_handle        = @(~,~) load_several_experiments(roi, data_folders, mask_name);
    uicontrol('Parent',hFig,'Style','pushbutton','String','Load ROI','Units','normalized','Position',[0.8 0.95 0.15 0.05],'callback',f_handle);
    
    br = get_id(data_folders{1},0,roi);
    rois = get_id(data_folders{1},br,0);
    f_handle2        = @(~,~) load_several_experiments(rois, data_folders, mask_name);    
    uicontrol('Parent',hFig,'Style','pushbutton','String','Load Branch','Units','normalized','Position',[0.8 0.90 0.15 0.05],'callback',f_handle2);
end