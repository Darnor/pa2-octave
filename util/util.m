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

function v = to_vector(A, l, b)
  for i = 1:b
    for j = 1:l
      v((i-1)*l+j, 1) = A(i,j);
    endfor
  endfor
endfunction

function A = to_matrix(v, l, b)
  for i = 1:b
    for j = 1:l
      A(i,j) = v((i-1)*l+j);
    endfor
  endfor
endfunction
