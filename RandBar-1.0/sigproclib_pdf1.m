
% -----------------------------------------------------------------
%  sigproclib_pdf1.m
%
%  This functions computes the probability density funtion
%  of a given data.
%
%  input:
%  data     - data vector
%  numbins  - number of bins
%  data_min - minimum value for bins interval (optional)
%  data_max - maximum value for bins interval (optional)
%
%  output:
%  bins  - bins locations vector
%  freq  - frequency counts vector
%  area  - area under the histogram (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [bins,freq,area] = sigproclib_pdf1(data,numbins,data_min,data_max)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    elseif nargin == 3
        error('Enter a value for xmax.')
    elseif nargin == 2
    	data_max = max(data);
    	data_min = min(data);
    end
    
    % check arguments
    if ( data_max < data_min )
        error('xmax must be greather than xmin.')
    end
    
    
	% compute data vector size
	Ndata = length(data);

	% compute bin width
	binwidth = (data_max-data_min)/(numbins-1);

	% compute bins locations vector
	bins = (data_min:binwidth:data_max)';

	% compute frequency counts vector
	freq = histc(data,bins);

	% normalize frequency counts vector
	freq = freq./(Ndata*binwidth);

	% compute area under the histogram
	area = binwidth*sum(freq);

return
% -----------------------------------------------------------------
