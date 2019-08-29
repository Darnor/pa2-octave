%%
% function functions
function delta = dirac_delta(x)
  for k=1:size(x,2)
    if (x(k) < 1e-16 && x(k) > -1e-16)
      delta(k) = 1;
    else
      delta(k) = 0;
    endif
  endfor
  delta = delta';
endfunction

function r = periodic_fun(N)
  x = linspace(0, 2*pi, N);
  r = 10 + sin(x) + 2 * cos(10 * x);
  r = r';
endfunction

