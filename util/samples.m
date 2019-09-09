%%
% params: x
%
% Dirac Delta
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

%%
% params: N
function r = sample_periodic_fun(N)
  x = linspace(0, 2*pi, N);
  r = 10 + sin(x) + 2 * cos(10 * x);
  r = r';
endfunction

%%
% params: im, J = 10, cutoff = 20
function [T, ind] = sample_image_wavelets(im, J = 10, cutoff = 20)
  ind = get_local_maxima_indexes(im, @(x) cutoff_band_filter(x, cutoff));
  L = laplace_from_indexes(ind, @image_pixel_euclidean_distance);
  [chi, lambda] = eig(L);
  lambda = diag(lambda);
  T = wavelets(@g, @h, chi, lambda, J);
endfunction

%%
% params: C, T, ind
function w = sample_image_W(C, T, ind)
  for i = 1:size(C, 2)
    f = function_from_indexes(C{i}, ind);
    w{i} = W(f, T);
  endfor
endfunction

%%
% params: w, j
function sample_plot_W(w, j)
  y_max = max(w{1}{j});
  y_min = min(w{1}{j});
  for k = 2:size(w, 2)
    y_max = max(max(w{k}{j}), y_max);
    y_min = min(min(w{k}{j}), y_min);
  endfor
  for k = 1:size(w, 2)
    figure(k);
    plot(w{k}{j});
    axis ([0 size(w{k}{j}, 1)-1 y_min y_max]);
  endfor
endfunction
