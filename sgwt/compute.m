source("util/util.m")

%%
% params: C, l, b, x0, y0, T
%
% Compute image wavelet difference using a lxb window
function d = compute_wavelet_diff(C1, C2, l, b, x0, y0, T)
  Q1 = C1(y0+1:y0+b,x0+1:x0+l);
  Q2 = C2(y0+1:y0+b,x0+1:x0+l);

  f1 = to_vector(Q1, l, b);
  f2 = to_vector(Q2, l, b);

  meansq(f1-f2)

  w1 = W(f1, T);
  w2 = W(f2, T);

  for i = 1:size(w1, 2)
    d{i} = w1{i}-w2{i};
  endfor
endfunction
