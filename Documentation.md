# Analyse a recording (arboreal_scan objects)

TODO

This will be the doc for the Arboreal_Scan objects



## Build the arboreal_scan object

```matlab
data_folder = 'D:\Curated Data\2019-10-03\experiment_5\20-36-43' 

%% Create the arboreal_scan object for data foldr
a_s = arboreal_scan(data_folder);

%% Prepare tree information. 
% Option 1, settings.txt is in the TOP_FOLDER
a_s.prepare_tree();
% Prepare tree information. Option 2, provide the path to settings.txt
a_s.prepare_tree('path/to/settings.txt');
% Prepare tree information. Option 2, provide the batch_params structure
a_s.prepare_tree(batch_params);

%% Prepare the extraction (concatenation etc...). This uses default values (if you type a_s.prepare_extraction() the default values are printed in an error message)
a_s.prepare_extraction(data_folder);

%% Save
a_s.save();
```



##  Basic Use

```matlab
%% Plot ROIs (for advanced use of this method, see dedicated doc)
a_s.plot_value_tree

%% Plot soma population 
a_s.plot_population();

%% Plot the distance from soma (use to check that branches are well reconnected)
a_s.plot_dist_tree
```





# Analyse an experiment (arboreal_scan_experiment objects)

Load an `arboreal_scan_experiment` object, that you have previously generated. This object regroups multiple `arboreal_scan` objects that shared the same tree structure. The loaded object is called "obj" by default.

[TODO : Link to how they were generated ]

![image-20210604171843850](media/Documentation/image-20210604171843850.png)

## Processing pipeline

Essentially, the analysis process consist in :

- Loading all the relevant arboreal_scan objects that come from the same cell into an arboreal_scan_experiment
- Rescaling traces so you can compare dim and bright regions
- Identify Global/Local events
- Flag problematic/low quality/wrongly selected ROIs so they will be ignored
- Regroup traces into "bins", i.e. regions that will be processed together (e.g. basolateral vs apical, or bins at some set distance from the soma)
- Process data by bin or by ROI. Among these processing options you can
  - Analyse correlation or similarities between regions
  - Extract events
  - Do some dimensionality reduction analysis on each region or ROI
  - Correlate against behaviours

#### Loading/object building 

This function is typically controlled by the constructor, although you can update the object later if you changed the list of recordings, or if you chaged the extraction method (e.g. new post-hoc MC, different mask...)

```matlab
%% Case 1
%% To build the object from extracted arboreal_scans (recommanded)
obj = arboreal_scan_experiment(source_folder); % where source_folder is a folder containing multiple extracted arboreal_scans 												   % objects (files can be in subfolders) 
%% Case 2
	%% To build the object from data_folders
obj = arboreal_scan_experiment(source_folder); % In the absence of extracted arboreal_scan, The code will ask you if you want to process the data directly. Note that as for regular arboreal_scan extraction, the folder must follow the standard day/expe/data_folder structure, and there must be a settings.txt file indicating how to reconnect the tree. see arboreal_scan documentation for details

	%% To build the object from data_folders using specific analysis_params
obj = arboreal_scan_experiment(source_folder, analysis_params('smoothing', 20));

	%% To build the object from data_folders using a custom settings.txt file
obj = arboreal_scan_experiment(source_folder, '', 'path/to/settings.txt'');

%% Case 3
%% To update the object if you added/removed arboreals scans
obj.update(); % Note that since this will affect the binned traces, all analysis need to be regenerated. All analysis fields will be cleared
```

#### Analysis can be handled by `obj.process()`

```matlab
%% Process the entire tree at once
obj.process();

%% Process using one of the supported binning method
obj.process({'depth',100}); % by depth, using bins of 100 um

%% Process without generating figures
obj.rendering = false;
obj.process(...);
```

#### Save the result

```matlab
%% To save the object
obj.save()

%% If you want to auto-save results, you can do the following
obj.auto_save_analysis = true;
obj.process(...);

%% To save the figures
obj.save_figures();

%% or, to auto-save the figures, you can do the following
obj.auto_save_figures = true;
obj.process(...);
```

## Object structure

```matlab
            source_folder: 'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER'
     extracted_data_paths: {1×3 cell}
           arboreal_scans: {[1×1 arboreal_scan]  [1×1 arboreal_scan]  [1×1 arboreal_scan]}
                      ref: [1×1 arboreal_scan]
                timescale: [1×1 struct]
        global_median_raw: [1485×1 single]
     global_median_scaled: []
                        t: [1×1485 double]
             batch_params: [1×1 struct]
              need_update: [1 1 1]
        extraction_method: 'median'
         extracted_traces: {[495×107 single]  [495×107 single]  [495×107 single]}
    extracted_traces_conc: [1485×107 single]
                   n_ROIs: 107
                     demo: 0
               filter_win: [5 0]
              filter_type: 'gaussian'
                 peak_thr: 2
                  cc_mode: 'peaks'
          rescaled_traces: [1485×107 single]
                  detrend: 1
           rescaling_info: [1×1 struct]
              binned_data: [1×1 struct]
                    event: [1×1 struct]
            event_fitting: [1×1 struct]
              variability: [1×1 struct]
                crosscorr: 1
           dimensionality: [1×1 struct]
       external_variables: {[1×1 struct]  [1×1 struct]  [1×1 struct]}
              spiketrains: []
               behaviours: [1×1 struct]
           default_handle: @(x)load_several_experiments(x,cellfun(@(x)x.data_folder,obj.arboreal_scans,'UniformOutput',false),use_mask)
             bad_ROI_list: [1 2 3 4 5 20 27 28 31 32 33 34 47 50 72 106 107]
                rendering: 1
```

## Some general properties

`obj.source_folder` indicates where these individual arboreal_scans recordings were located, and `extracted_data_paths` point at the exact file used for loading these data.

`obj.arboreal_scans` is a cells array containing the individual arboreal_scans you extracted, i.e. an object containing information about the tree structure, the experiment (header, behavioural data) and the extracted calcium signals, using specific masks and registration options. 

> Note : if you were to change the registration method or the masks, you will have to re-extract the calcium signals, and therefore regenerate the arboreal_scans objects and any object/analysis that depend on it (e.g.  the arboreal_scan_experiment object). The extraction options initially used can be obtained in `obj.arboreal_scans{1}.analysis_params`. The mask used can be seen in `obj.arboreal_scans{1}.mask`

> Note : obj.arboreal_scans{N}  DOES NOT contain as much data than the source arboreal_scan object indicated (i.e. one value per voxel along along a ROI), but ONLY one trace per ROI. This value is obtained through a specific compression mechanism. For more details about the compression mechanism, see the [corresponding section](#....). 

`obj.ref` points at the first arboreal_scan, to simplify function calls. This assume that all your arboreal_scans are identical

`obj.timescale` return, for each recording 

* the imaging sampling rate (`obj.timescale.sr`) ; measured if possible, estimated otherwise.
* The number of imaging timepoints `obj.timescale.tp`
* The duration (`obj.timescale.durations`) ; measured if possible, estimated otherwise.
* A timescale array for each recording, starting at 0 for each recording `obj.timescale.rec_timescale`
* A global, continuous timescale (no gaps between trials or recordings), where all records are contactenated `obj.timescale.global_timescale`
* The corresponding t_start for each one of these recordings `obj.timescale.t_start_nogap`
* *TODO: THE EXACT, REAL T START*

`obj.t` is a pointer to `obj.timescale.global_timescale`



`obj.batch_params` indicates the arboreal_scan structure info for the cells that are required to rebuild the tree correctly with the Trees toolbox (soma depth, primary branches, manual reconnections, excluded branches, pia and soma depth). This is actually a pointer from `obj.arboreal_scans{1}.batch_params`, and they were obtained at the time of the arboreal_scans extraction. See the main documentation on that topic.

`obj.need_update` .. TODO

`obj.n_ROIs` : the number of segments per tree, and therefore the number of traces

`obj.extracted_traces` : the extracted traces originating from each arboreal_scan, using the compression method defined in [XXXXXX]. This a [1 x N_records] cell array, of {T x N_ROI} cells. `obj.extracted_traces_conc` corresponds to the same traces, concatenated in time.

## Rescaling

`obj.rescaled_traces`. This is obtained from `obj.extracted_traces`, but where each ROIs were rescaled in regard to each other, so that the median amplitude of the events in each ROI is identical to the cells median. The rescaling is done by calling `obj.rescale_traces()`, which generate the `obj.rescaling_info` field. Each ROI is rescaled using a unique offset and scalar across all the recordings.

Rescaling figures 1012 (Scaling Factor per ROI, per trial) and 1013 (Global scaling factor and offset per ROI) can be obtained by calling `obj.plot_rescaling_info`

<img src="media/Documentation/image-20210604180713086.png" alt="image-20210604180713086" style="zoom: 33%;" />

<img src="media/Documentation/image-20210604180748594.png" alt="image-20210604180748594" style="zoom: 33%;" />



## Binning

In order to analyse signal differences across different parts of the dendritic tree, you can regroup ROIs based on multiple morphometric criteria, such as the depth, the distance from the soma etc...

Most analyses are done on a group-by-group basis, so one of the first steps after rescaling the traces is to choose the relevant binning method. To generate the groups, call `obj.prepare_binning(CONDITION)`, where the condition can be one of the preset modes defined in the `arboreal_scan.get_ROI_groups()` method, which uses `get_ROIs_subsets()`, or a custom segmentation.

To pass a preset mode, the CONDITION must be a STR 'condition', or a CELL ARRAY {'condition', metric} matching one of the following : 

* {'distance',DISTANCE_BIN_SIZE_UM}
* {'depth',DEPTH_BIN_SIZE_UM}. Note that bin values are centred around the soma, not from the pia, so that basolateral and apical dendrites are usually segregated.
* {'random', N_ROI_PER_BIN}
* 'branch' ; Each branch (as they were traced) form a group
* 'order' ; Branch order (primary, secondary etc..)

For a custom control, you can pass an array of cells, each cell will form a group and the condition will be labelled as 'custom'.

```matlab
%% Example 1 : Segment by bits of 100um from the soma
obj.prepare_binning({'depth',100})

%% Example 2 : Use customized groups
obj.prepare_binning([{1:10},{11:20},{21:30}]);

%% Example 3 : No binning (all ROIs in 1 group)
obj.prepare_binning();
```

The condition is stored in `obj.binned_data.condition`, the bin values are stored in `obj.binned_data.metrics`, the ROIs for each group are in `obj.binned_data.groups` (same order) and the legends displayed on figures are in `obj.binned_data.groups`

## Computing Median Traces

`obj.binned_data.global_median` correspond the median trace of all recordings, and somewhat constitute a signal representative of the entire tree. The  median trace is obtained by calling `median_trace = obj.set_median_traces()`. By default, this uses the median of the rescaled traces, although you can use the median of the raw traces instead (`median_trace = obj.set_median_traces(true)`).

`obj.binned_data.median_traces` correspond to the median trace PER GROUP, as defined by your binning method (see [binning doc](#Binning) for `obj.prepare_binning()`).

> Median traces an be displayed without or with gaussian smoothing

```matlab
obj.plot_median_traces; 			% basic median
obj.plot_median_traces(20); 		% median with a 20 point symetrical smoothing
obj.plot_median_traces([20, 0]);    % median with a 20 point asymetrical smoothing
```

## Event Detection

You can detect large transients, and store they time of occurrence in `obj.event`. 

This uses the `obj.detect_events()`. Event detection can be using either a peak_amplitude approach [TODO : TO PUT BACK], or an approach that looks at correlated signal variations across the tree (either all of it, or a selected region such as the peri somatic area. see `idx_filter` input).

The default uses the following steps:

- Pairwise correlations between each pair of ROI are computed, using a moving correlation window (see `corr_window` input). If no window_side is given, the size is set as the average event width (as found by the native matlab `findpeaks` function, which returns here the width at half prominence of all events that are 2x signal RMS).

- The average of all these Pairwise correlations is created, and indicate how correlated is the signal across the tree. A value of 1 indicates that an event covaried across every ROI of the tree, while a value of 0 suggest completely random variations (e.g. no activity). As some regions may be belonging to other neurons, the maximum may not reach 1. Therefore, the result is renormalized to the maximal values

  > These uncorrelated ROIs can be detected and excluded if we consider that they are uncorrelated because belonging to another cell, or are too faint to provide a usable signal. see XXXX

- All events that are at least present in 10% of the tree are selected as "events" This threshold can be modified by using the `cutoff` input. 

  - The time of each detected event is stored in `obj.event.t_corr`. Note that this is using the correlation window, and the actual time of peak may vary and is stored in `obj.event.t_peak` (see below).
  - The amplitude of the peak indicate how tree-wide is the correlation, and is stored in `obj.event.globality_index`
  - "Global" Events are arbitrarily set as events propagating in more than 50% of the tree, and are indicated in `obj.event.is_global`
  - Original "raw" amplitude of the fluorescent signal is stored in `obj.event.peak_v`
  - The lower and upper time range of each event are stored in `obj.event.t_win` (as defined by the onset from/offset to baseline, or by the end/start of another surrounding event) 

If you want to detect somatic events only, you can use `obj.detect_events('soma')` which will only use the ROIs located the closest to the soma (on primary branches). If there is no soma (e.g. L5 cell, then the most proximal segments will be used instead)



ADD FIGURES



### Signal Compression

`arboreal_scan_experiments` have one trace per ROI. This means that each *"2D"* X-Y-T patch is compressed into a single *"0D"* 1-1-T signal array. The compression along the patch axis (that converts a X-Y-T patch into a *"1D"* X-1-T line) happens during the initial signal extraction, and these values are stored in individual `arboreal_scan` objects. This step can not be modified without re-extracting the arboreal_scans. The compression that converts *"1D"* X-1-T lines into *"0D"* 1-1-T array happens when first generating the arboreal_scan_experiment object, through the `obj.update_all_signals()` method (which actually call the `arboreal_scan.update_segment_signals()` method). The default approach takes the median signal, although this can be changed by passing a method parameter `obj.update_all_signals(method)`, with any method in {'min, 'max', median', 'mean'}. The current method is stored in the`obj.extraction_method` field.

> Reloading the data to change the compression method requires to point at the folder containing the `arboreal_scans`. This is done using `obj.extracted_data_paths`. If you moved you data, you'll need to update these paths

> Changing the method requires to reload all the individual `arboreal_scans`. To validate your changes, you need to save the object using `obj.save();`.

ADD IMAGE COMPRESSION

