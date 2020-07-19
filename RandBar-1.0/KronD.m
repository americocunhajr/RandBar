
% -----------------------------------------------------------------
%  KronD.m
%
%  This functions computes the Kronecker delta tensor.
%
% input:
% i        - 1st index
% j        - 2nd index
%
% output:
% delta_ij - delta_ij component
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 30, 2011
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function delta_ij = KronD(i,j)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    
    if i == j
        delta_ij = 1;
    else
        delta_ij = 0;
    end

return
% -----------------------------------------------------------------
