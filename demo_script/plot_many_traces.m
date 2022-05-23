%% TODO: ADD DOC, SET SPACING FACTOR AS A INPUT PARAMETER

function [traces, fig, offsets] = plot_many_traces(traces, fig, varargin)
    if nargin < 2  || isempty(fig)
        fig         = figure(); hold on;
    elseif isnumeric(fig) || ishandle(fig)
        if isa(fig, 'matlab.graphics.axis.Axes') % subplot
            subplot(fig)
        else % fig
            fig     = figure(fig); hold on;
        end
    else
        fig = [];
    end
    %traces = smoothdata(squeeze(traces), 'gaussian',[5, 0]);
    
    if ndims(traces) > 2
        traces      = squeeze(traces);
    end
    
    spacing_factor  = 1    
    spacing         = nanmedian(rms(traces))/spacing_factor;    
    m               = nanmin(traces(:));
    traces(:,1)     = m;
    traces(:,end)   = m;
    
    n_traces        = size(traces, 2);
    offsets         = linspace(0,n_traces*spacing,n_traces);    
    traces          = traces + offsets;
    if ~isempty(fig)
    	plot(1:size(traces, 1), traces, varargin{:});ylim([nanmin(traces(:)),nanmax(traces(:))])
    end
end

