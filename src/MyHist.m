function varargout = MyHist( data,binsize,varName,units)
%MYHIST Plot a histogram and adds mean/sd and median to the plot.
first  = floor(min(data/binsize)) * binsize;
last   = ceil(max(data/binsize)) * binsize;
[n,c] = hist(data, first : binsize : last);

if nargout < 2
    % If no or one output is requested, plot the histogram as a bar plot.
    handle = bar(c,n);
    if nargout == 1
        % Return the handle if one output argument is requested.
        varargout{1} = handle;
    end
    % Add title
    title(sprintf('%s, binsize = %.2f',varName,binsize))
    xlabel([varName ' (' units ')'])
    ylabel('count')
    text(0.95,0.95,{sprintf('mean (sd) = %.2f (%.2f)',nanmean(data),nanstd(data)),...
        sprintf('median = %.2f%.2f)',nanmedian(data)),...
        sprintf('n = %d',sum(~isnan(data)))},...
        'HorizontalAlignment','Right','VerticalAlignment','Top',...
        'Color','k','Units','Normalized')
elseif nargout == 2
    % If two outputs are requested, don't plot the data but return the bin
    % centres (c) and counts per bin (n) as outputs.
    varargout{1} = c;
    varargout{2} = n;
    
end

end % of function

