
% -----------------------------------------------------------------
%  graph_semilogy4.m
%
%  This functions plots semi log y graph with four curves.
%
%  input:
%  x1     - x data vector 1
%  y1     - y data vector 1
%  x2     - x data vector 2
%  y2     - y data vector 2
%  x3     - x data vector 3
%  y3     - y data vector 3
%  x4     - x data vector 4
%  y4     - y data vector 4
%  gtitle - graph title
%  leg1   - legend 1
%  leg2   - legend 2
%  leg3   - legend 3
%  leg4   - legend 4
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
%  last update: Fev 20, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_semilogy4(x1,y1,x2,y2,x3,y3,x4,y4,...
                               gtitle,leg1,leg2,leg3,leg4,...
                               xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)
    
    % check number of arguments
    if nargin < 20
        error('Too few inputs.')
    elseif nargin > 21
        error('Too many inputs.')
    elseif nargin == 20
        flag = 'none';
    end

    % check arguments
    if length(x1) ~= length(y1)
        error('x1 and y1 vectors must be same length')
    end
    
    if length(x2) ~= length(y2)
        error('x2 and y2 vectors must be same length')
    end
    
    if length(x3) ~= length(y3)
        error('x3 and y3 vectors must be same length')
    end
    
    if length(x4) ~= length(y4)
        error('x4 and y4 vectors must be same length')
    end
    
    if length(x1) ~= length(x2)
        error('x1 and x2 vectors must be same length')
    end
    
    if length(x1) ~= length(x3)
        error('x1 and x3 vectors must be same length')
    end    
    
    if length(x1) ~= length(x4)
        error('x1 and x4 vectors must be same length')
    end
    
    if length(x2) ~= length(x3)
        error('x1 and x3 vectors must be same length')
    end
    
    if length(x2) ~= length(x4)
        error('x2 and x4 vectors must be same length')
    end
    
    if length(x3) ~= length(x4)
        error('x3 and x4 vectors must be same length')
    end
    
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = semilogy(x1,y1,'-b');
    hold all
    fh2 = semilogy(x2,y2,'-.r');
    fh3 = semilogy(x3,y3,':m');
    fh4 = semilogy(x4,y4,'--k');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',20);
    title(gtitle,'FontSize',20,'FontName','Helvetica');
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
    
    set(fh1,'LineWidth',2.0);
    set(fh1,'MarkerSize',4.0);
    set(fh1,'MarkerFaceColor','w');
    set(fh1,'MarkerEdgeColor','b');
    set(fh2,'LineWidth',2.0);
    set(fh2,'MarkerSize',4.0);
    set(fh2,'MarkerFaceColor','w');
    set(fh2,'MarkerEdgeColor','r');
    set(fh3,'LineWidth',2.0);
    set(fh3,'MarkerSize',2.0);
    set(fh3,'MarkerFaceColor','w');
    set(fh3,'MarkerEdgeColor','k');
    set(fh4,'LineWidth',2.0);
    set(fh3,'MarkerSize',2.0);
    set(fh3,'MarkerFaceColor','w');
    set(fh4,'MarkerEdgeColor','k');
    leg = legend(leg1,leg2,leg3,leg4,3);
    set(leg,'FontSize',12);
    set(leg,'interpreter', 'latex');
    labX = xlabel(xlab,'FontSize',20,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',20,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------
