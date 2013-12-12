function rc = f_rate(vecnum, id)

% f_rate.m
% Estimates the convergence rate of a given sequence
% of scalars (vecnum=Vector) or vectors (vecnum=Matrix).
%
% Usage: K = f_rate ( vecnum, id)
% ==or== K = f_rate ( vecnum )
%
% where vecnum is of size (m x n), with n = dimension, and
% id=1 is to show all intermediate estimates.
% For usage 2, id=0;
% in all cases the very last estimate is shown.
% For example, to check the convergence of the iterates,
% (1, 2, 3, 4, 4.5, 5, 5.11)
% USE vecnum = [1 2 3 4 4.5 5 5.11]';
% kc_1 = f_rate(vecnum,1) %% Giving 1.3914 ... 2.8614
% kc_2 = f_rate(vecnum) %% Giving 1.7101 and 2.8614

if (nargin < 1), help f_rate; return, end
if nargin==1, id=0; end
m = size(vecnum,1);
abort_id = 0;
x_exact = vecnum(m,:) ; % Take as exact solution "x"
keep = [] ; rc = []; i=6;
while (i ~= m)
    d_n_m = norm(vecnum(i-2,:)-x_exact) ; % d_(n-1)
    d_n = norm(vecnum(i-1,:)-x_exact) ; % d_n
    d_n_p = norm(vecnum(i,:) -x_exact) ; % d_(n+1)
    if (d_n_m < eps)
        d_n_m = eps*7;
    end
    if (d_n < eps)
        d_n = eps*7;
    end
    if (d_n_m < eps)
        d_n_m = eps*7;
    end
    if (d_n < eps)
        d_n = eps*7;
    end
    if (d_n_p < eps)
        d_n_p = eps*7;
        i = m - 1; 
        abort_id = 1; % No need to continue
    end
    top = d_n_p / d_n ; % numerator
    bottom = d_n / d_n_m ; % denominator
    if (top == 1)
        top = 1 + eps*7 ;
    end % Safeguard for Log
    if (bottom == 1)
        bottom = 1 + eps*7 ;
    end
    rc = log(top)/log(bottom) ; % Estimate of the order
    if (abort_id == 0 && rc > 0), keep = [keep; rc]; end
    i = i + 1;
end % End WHILE-loop
disp(' ');

if id == 1
    rc = keep ;
    disp(' (Individual estimates for rate)') ;
else
    m = length(keep); m=keep(m);
    fprintf('Average Rate =%5.1f\n',mean(keep));
    fprintf('Last estimate=%5.1f\n',m);
    rc = [mean(keep) m];
end % End for IF