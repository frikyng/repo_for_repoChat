function save_myfig(this_fig,file_name,file_ext)
% Saves <this_fig> to <file_name> as a <file_ext>
%   If multiple extensions are given, saves as each format
%   If extension is not recognized, it tries anyway
%   'pdf_small' creates a low-resolution PDF
% make sure the figure has updated
drawnow
ext_options = {'fig','pdf','png','jpg','jpeg','eps','pdf_small'};
%% Parse inputs
% if no handle given, save the current figure
if nargin<1 || isempty(this_fig)
    this_fig=gcf;
end
% if no name given, call it "figure.<ext>"
if (nargin<2 || isempty(file_name))
    file_name='figure';
% trim off existing extensions
elseif ~isempty(file_name) && endsWith(file_name,ext_options)
    [filepath,name,~] = fileparts(file_name);
    file_name = [filepath filesep name];
end
% make sure the arguments are in the right format
file_name = convertStringsToChars(file_name);
if nargin<3 || isempty(file_ext)
    file_ext={};
elseif ischar(file_ext)
    file_ext={file_ext};
elseif isnumeric(file_ext)
    file_ext=ext_options(file_ext);
end
%% check for each known extenion and save in that format
% JPG
if any(strcmpi(file_ext,'jpg')) || any(strcmpi(file_ext,'jpeg'))
    print(this_fig,'-djpeg',[file_name '.jpg'])
end
% PDF (low resolution)
if any(strcmpi(file_ext,'pdf_small'))
    % I always print in landscape, but this may not make sense with your figures
    orient(this_fig,'landscape');
    print(this_fig,'-dpdf','-r72','-bestfit',[file_name '.pdf']);
    orient(this_fig,'portrait');
end
% PNG
if any(strcmpi(file_ext,'png'))
    print(this_fig,'-dpng',[file_name '.png'])
end
% EPS
if any(strcmpi(file_ext,'eps'))
    print(this_fig,'-depsc',[file_name '.eps'])
end
% PDF (normal resolution)
if any(strcmpi(file_ext,'pdf'))
    % I always print in landscape, but this may not make sense with your figures
    orient(this_fig,'landscape');
    tmp=get(this_fig,'Position');
    set(this_fig,'PaperUnits','points','PaperSize',tmp(3:4))
    print(this_fig,'-dpdf','-r300',[file_name '.pdf']);
    orient(this_fig,'portrait');
end
if any(strcmpi(file_ext,'fig'))
    % clear warnings for the try-catch below
    lastwarn('')
    % check whether the figure is visible
    invisible = false;
    if ~strcmp(get(this_fig,'visible'),'on')
        invisible = true;
        set(this_fig,"Visible","on");
    end
    % try the normal way
    try
        saveas(this_fig, [file_name '.fig'],'fig')
    catch
        % warning('Variable ''d'', is larger than 10GB and could not be saved')
    end
    warnMsg = lastwarn;
    % use v7.3 saving for large data
    if ~isempty(warnMsg)
        % warning('Fig Variable too large for v7, using v7.3 saving')
        hgsave(this_fig, [file_name '.fig'], '-v7.3');
        fprintf(['Fig actually saved\n'])
    end
    if invisible
        set(this_fig,"Visible","off");
        invisible = false;
    end
end
% Try unrecognized file extensions
for ind = 1:length(file_ext)
    if ~any(strcmp(file_ext{ind},ext_options))
        print(this_fig,['-d' file_ext{ind}],[file_name '.' file_ext{ind}])
    end
end
end