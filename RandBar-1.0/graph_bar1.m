
% -----------------------------------------------------------------
%  graph_bar1.m
%
%  This functions plots a bar histogram with one series.
%
%  input:
%  bins1  - bins locations vector
%  freq1  - frequency counts vector
%  gtitle - graph title
%  xlab   - x axis label
%  ylab   - y axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Mar 21, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_bar1(bins1,freq1,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)
	
    % check number of arguments
    if nargin < 10
        error('Too few inputs.')
    elseif nargin > 11
        error('Too many inputs.')
    elseif nargin == 10
        flag = 'none';
    end

    % check arguments
    if length(bins1) ~= length(freq1)
        error('bins1 and freq1 vectors must be same length')
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = bar(bins1,freq1,0.8);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
    if ( strcmp(xmin,'auto') || strcmp(xmax,'auto') )
        xlim('auto');
    else
        xlim([xmin xmax]);
    end

    if ( strcmp(ymin,'auto') || strcmp(ymax,'auto') )
        ylim('auto');
    else
        ylim([ymin ymax]);
    end
    
    set(fh1,'FaceColor','none');
    set(fh1,'EdgeColor','b');
    set(fh1,'LineStyle','-');
    labX = xlabel(xlab,'FontSize',16,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',16,'FontName','Helvetica');
    %set(Xlab,'interpreter','latex');
    %set(Ylab,'interpreter','latex');
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------
