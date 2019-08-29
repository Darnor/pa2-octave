%%
% params: A, value
%
% Check if all values in A match the value
function r = all_match(A, value)
  r = 1 && all(A == value);
endfunction

%%
% params: A, value
%
% Check if all values in A match the value
function r = any_match(A, value)
  r = any(any(A == value) == 1);
endfunction
