
%% Fixes and standardizes the names of labels
%  The function standardizes the naming convention of labels in a dataset,
%  especially useful when labels may have variations due to the use of
%  different naming protocols.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[labels] = fix_labels(labels)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	labels(Array of Strings OR Categorical Array):
%                                   The labels to be standardized. The
%                                   labels can be a cell array of strings
%                                   or a categorical array.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	labels(Array of Strings OR Categorical Array)
%                                   The labels after the naming has been
%                                   standardized.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * Labels are renamed according to a fixed pairs list defined in the
%   function. If the label is not found in the pairs list, it is kept as is.
% -------------------------------------------------------------------------
%% Examples:
% * Standardize a cell array of string labels
% 	labels = fix_labels({'encoder', 'EyeCam_L_forelimb'});
%
% * Standardize a categorical array of labels
% 	labels = fix_labels(categorical({'encoder', 'EyeCam_L_forelimb'}));
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% Paste license here.
% -------------------------------------------------------------------------
% Revision Date:
% 	09-06-2023
% -------------------------------------------------------------------------
% See also: 

function labels = fix_labels(labels)
    if iscategorical(labels)
        oldnames = categories(labels);
    else
        oldnames = labels;
    end
    
    if any(contains(oldnames,'\_'))
        pairs = {   'encoder', 'Running Speed';...
                    'EyeCam\_L\_forelimb', 'I. Forelimb MI' ; ...
                    'EyeCam\_R\_forelimb', 'C. Forelimb MI' ; ...
                    'BodyCam\_L\_whisker', 'I. Whisker MI' ; ...
                    'BodyCam\_R\_whisker', 'C. Whisk.Pad MI' ; ...
                    'EyeCam\_Perioral', 'Peri-oral MI' ; ...
                    'trigger', 'C. AirPuff' ; ...
                    'baseline', 'F0'};
    else
        pairs = {   'encoder', 'Running Speed';...
                    'EyeCam_L_forelimb', 'I. Forelimb MI' ; ...
                    'EyeCam_R_forelimb', 'C. Forelimb MI' ; ...
                    'BodyCam_L_whisker', 'I. Whisker MI' ; ...
                    'BodyCam_R_whisker', 'C. Whisk.Pad MI' ; ...
                    'EyeCam_Perioral', 'Peri-oral MI' ; ...
                    'trigger', 'C. AirPuff' ; ...
                    'baseline', 'F0'};
    end

    newnames = oldnames;
    for el = 1:numel(oldnames)
        if ~contains(oldnames{el}, 'shuffle')
            newnames{el} = pairs{find(~cellfun(@isempty, (strfind(pairs, oldnames{el})))), 2};
        end
    end

    if iscategorical(labels)
        labels          = renamecats(labels,oldnames,newnames);
    else
        labels          = newnames;
    end
end

