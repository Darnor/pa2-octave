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

function print_wavelets(w)
  for j = 1:size(w, 2)
    O = (w{j}{1}).^2;
    for k = 2:size(w{j}, 2)
      O = [O (w{j}{k}).^2];
    endfor
    O = O';
    name = ["~/abs_value_" num2str(j) ".svg"];
    imagesc(O);
    colorbar();
    print("-dsvg", name);
  endfor
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
function [T, ind, d] = sample_image_wavelets(im, indexFn, J = 10, max_N = 1000)
  ind = indexFn(im, max_N);
  d = image_pixel_euclidean_distance_all(ind, true);
  L = laplace_from_distances(d);
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
    plot(abs(w{k}{j}).^2);
    axis ([0 size(w{k}{j}, 1)-1 y_min y_max]);
  endfor
endfunction

%%
% params: orig_im, s = 1, n = 3, fps = 50, J = 10, max_N = 1000, fn = @(x, y) 20*sin(5*x-y)
function [S, phases, T, ind, d, w] = sample_first_n_stretched_images_fps(orig_im, s = 1, n = 3, fps = 50, J = 10, max_N = 1000, fn = @(x, y) 20*sin(10*x-y))
  if (n > fps)
    n = fps;
  endif
  phases = linspace(0, 2*pi, fps);
  if (n < 1)
    S{1} = orig_im;
  endif
  A = stretch_modulate_image(orig_im, @(x) fn(x, phases(1)));
  for k = 1:n-(s-1)
    S{k} = stretch_modulate_image(orig_im, @(x) fn(x, phases(k+(s-1))));
  endfor
  [T, ind, d] = sample_image_wavelets(A, J, max_N);
  w = sample_image_W(S, T, ind);
  N = size(w{1}, 2);
  for j = 1:N
    %sample_plot_W(w, j, (j-1)*size(w, 2));
  endfor
  offset = N*size(w, 2);
  for j = N+1:N+n
    %figure(offset + (j - N));
    %sample_mesh_image_with_indices(S{j-N}, ind);
  endfor
endfunction

%%
% run full image transform sample
function sample_full()
  C = read_all_images("images");
  [T, ind, d] = sample_image_wavelets(C{1});
  w = sample_image_W(C, T, ind);
  for j = 1:size(w{1}, 2)
    sample_plot_W(w, j, (j-1)*size(w, 2));
  endfor
endfunction

%%
% params: img, ind
function sample_mesh_image_with_indices(img, ind)
  G = index_to_graph(ind, size(img, 1), size(img, 2));
  imagesc(img+200*G);
endfunction

%%
% params: img, ind, d
function sample_mesh_image_with_graph(img, ind, d)
  sample_mesh_image_with_indices(img, ind);
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
