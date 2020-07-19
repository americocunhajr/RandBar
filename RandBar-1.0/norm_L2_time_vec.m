
% -----------------------------------------------------------------
%  norm_L2_time_vec.m
%  
%  This function computes the L2 norm of the 2-norm of the 
%  time dependent vector t -> Q(t) in R^{Ndofs}.
%  
%                -- t1
%               |
%  normL2_Q :=  | ||Q(t)||^2 dt
%               |
%             -- t = t0
%  
%  
%  Input:
%  t0 - initial time
%  t1 - final   time
%  Q  - (Ndofs x Ndt) time dependent vector
%
%  Output:
%  normL2_Q - L2 norm of the 2-norm of Q(t)
% -----------------------------------------------------------------
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jan 17, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function normL2_Q = norm_L2_time_vec(t0,t1,Q)
    
    % check number of arguments
    if nargin < 1
        error(' Too few inputs.')
    elseif nargin > 3
        error(' Too many inputs.')
    end
    
    % check input arguments
    if ( t1 < t0 )
        error('t1 must be greather than t0.')
    end
    
    % compute the time dependent 2-norm of Q(t)
    norm2_Q = sqrt(sum(Q.^2));
    
    % compute the L2 norm of the 2-norm of Q(t)
    normL2_Q = norm_L2(t0,t1,norm2_Q);
    
return
% -----------------------------------------------------------------
