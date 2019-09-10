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
% params: im, J = 10, cutoff_max = 10, cutoff_min = -10
function [T, ind] = sample_image_wavelets(im, J = 10, cutoff_max = 10, cutoff_min = -10)
  ind_a = get_local_maxima_indexes(cutoff_min_filter(im, 10)); %cutoff_band_filter(x, cutoff_max, cutoff_min));
  ind_b = get_local_minima_indexes(cutoff_max_filter(im, -10));
  ind = [ind_a; ind_b];
  L = laplace_from_distances(image_pixel_euclidean_distance_all(ind, true, @(d) false));
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
% params: w, j, figure_id_offset
function sample_plot_W(w, j, figure_id_offset = 0)
  y_max = max(w{1}{j});
  y_min = min(w{1}{j});
  for k = 2:size(w, 2)
    y_max = max(max(w{k}{j}), y_max);
    y_min = min(min(w{k}{j}), y_min);
  endfor
  for k = 1:size(w, 2)
    figure(k+figure_id_offset);
    plot(w{k}{j});
    axis ([0 size(w{k}{j}, 1)-1 y_min y_max]);
  endfor
endfunction

%%
% run full image transform sample
function sample_full()
  C = read_all_images("images");
  D = C{2} - C{1};
  [T, ind] = sample_image_wavelets(D);
  w = sample_image_W(C, T, ind);
  for j = 1:size(w, 2)
    sample_plot_W(w, j, (j-1)*size(w, 2));
  endfor
endfunction

function sample_mesh_image_with_indices(img, ind)
  G = index_to_graph(ind, size(img, 1), size(img, 2));
  imagesc(img+200*G);
endfunction

%%
% params: img, ind, d
function sample_mesh_image_with_graph(img, ind, d)
  sample_mesh_image_and_graph(img, ind);
  hold on;
  for i = 1:size(ind, 1)
    for j = 1:size(ind, 1)
      if (d{i}(j) > 0)
        x = [ind(j, 2) ind(i, 2)];
        y = [ind(j, 1) ind(i, 1)];
        plot(x, y, 'r;;');
      endif
    endfor
  endfor
  hold off;
endfunction