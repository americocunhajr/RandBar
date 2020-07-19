
% -----------------------------------------------------------------
%  graph_contourf.m
%
%  This functions plots the contour map of a scalar function
%  F: R^2 -> R.
%
%  input:
%  x      - x mesh vector
%  y      - y mesh vector
%  F      - scalar field
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
function fig = graph_contourf(x,y,F,gtitle,xlab,ylab,...
                                    xmin,xmax,ymin,ymax,gname,flag)
                                    
    % check number of arguments
    if nargin < 11
        error('Too few inputs.')
    elseif nargin > 12
        error('Too many inputs.')
    elseif nargin == 11
        flag = 'none';
    end
    
    % generate 2D mesh grid
    xi      = linspace(xmin,xmax,length(x));
    yi      = linspace(ymin,ymax,length(y));
    [XI,YI] = meshgrid(xi,yi);
    FI      = griddata(x,y,F,XI,YI,'linear');

    fig = figure('Name',fname,'NumberTitle','off');
    
    fh = contourf(XI,YI,FI);
    colormap jet;
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
    labX = xlabel(xlab,'FontSize',16,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',16,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');

    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        graph_fixPSlinestyle(gname,gname);
    end


return
% -----------------------------------------------------------------