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

%%
% params: image, l, b, T
%
% Compute a windowed wavelet transform of the input image
function im = compute_windowed_wavelet_transform(image, l, b, T)
  [y, x] = size(image);

  N = y / b;
  M = x / l;

  for i = 1:N
    for j = 1:M
      ww = W(to_vector(image((i-1)*b+1:i*b, (j-1)*l+1:j*l), l, b), T);
      for k = 1:size(T, 2)
        R{k} = to_matrix(ww{k}, l, b);
      endfor
      Q{j} = R;
    endfor
    for k = 1:size(T, 2)
      for j = 1:M
        S{j} = Q{j}{k};
      endfor
      P{k} = cell2mat(S);
    endfor
    K{i} = P;
  endfor
  for k = 1:size(T, 2)
    for i = 1:N
      O{i} = K{i}{k};
    endfor
    im{k} = cell2mat(O');
  endfor
endfunction
