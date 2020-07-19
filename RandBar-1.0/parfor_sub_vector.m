
% -----------------------------------------------------------------
%  parfor_sub_vector.m
%
%  This function computes a sub vector from a given vector.
%
%  input:
%  v     - full vector
%  N1    - initial index
%  N2    - final index (optional)
%
%  output:
%  v_sub - sub vector
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 9, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function v_sub = parfor_sub_vector(v,N1,N2)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    elseif nargin == 2
        N2 = length(v);
    end
    
    % check arguments
    if N1 <= 0
        error('N1 must positive')
    end
    
    if N2 <= 0
        error('N2 must positive')
    end
    
    
    % compute sub vector
    v_sub = v(N1:N2);
    
return
% -----------------------------------------------------------------
