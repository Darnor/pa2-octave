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

function O = wavelet_image(w, ignore_h_kernel = true)
  O = (w{1}).^2;
  M = size(w, 2);
  if (ignore_h_kernel)
    --M;
  endif
  for k = 2:M
    O = [O (w{k}).^2];
  endfor
  O = O';
endfunction

function OO = wavelet_image_series(w, ignore_h_kernel = true)
  for j = 1:size(w, 2)
    OO{j} = wavelet_image(w{j});
  endfor
endfunction

function print_wavelet_image_compare(ww, w_num, type, ignore_h_kernel = true)
  y_max = 0;
  y_min = 0;
  l = 48;
  b = 27;
  OO = wavelet_image_series(ww, ignore_h_kernel);
  R = OO{1}(w_num, :);
  for j = 2:size(OO, 2)
    R = [R; OO{j}(w_num, :)];
  endfor
  y_max = max(R(:));
  y_min = min(R(:));
  print_wavelets(R, y_min, y_max, type, ["~/wavelet_compare_" num2str(w_num)]);
endfunction

function print_wavelets(O, y_min, y_max, type, file_name = "~/abs_value")
  imagesc(O, [y_min, y_max]);
  colorbar();
  print(["-d" type], [file_name "." type]);
endfunction

function print_wavelet_series(ww, type, ignore_h_kernel = true)
  y_max = 0;
  y_min = 0;
  OO = wavelet_image_series(ww, ignore_h_kernel);
  for j = 1:size(OO, 2)
    y_max = max(max(OO{j}(:), y_max));
    y_min = min(min(OO{j}(:), y_min));
  endfor
  for j = 1:size(OO, 2)
    print_wavelets(OO{j}, y_min, y_max, type, ["~/abs_value_" num2str(j)]);
  endfor
endfunction

function sample_trimmed_images(C, T, l, b, offset, ignore_h_kernel = true)
  ww = sample_image_series(C, T, @(im) to_vector(cut_image(im, l, b, offset), l, b));
  type = "png";
  print_wavelet_series(ww, type);
  N = size(ww{1}, 2);
  if (ignore_h_kernel)
    --N;
  endif
  for i = 1:N
    print_wavelet_image_compare(ww, i, type)
  endfor
  OO = wavelet_image_series(ww, ignore_h_kernel);
  % calculate the standard deviation of all wavelets for each image
  for i = 1:size(OO{1}, 1)
    R = OO{1}(i, :);
    for j = 2:size(OO, 2)
      R = [R; OO{j}(i, :)];
    endfor
    VV{i} = var(R')';
    SS{i} = std(R')';
    MM{i} = mean(R')';
  endfor

  V = cell2mat(VV);
  imagesc(V);
  colorbar();
  print(["-d" type], ["~/variance." type]);
  disp("VAR(VAR):");
  disp(var(V(:)));

  S = cell2mat(SS);
  imagesc(S);
  colorbar();
  print(["-d" type], ["~/stdeviation." type]);
  disp("STD(STD):");
  disp(std(S(:)));
  
  M = cell2mat(MM);
  imagesc(M);
  colorbar();
  print(["-d" type], ["~/mean." type]);
  disp("MEAN(MEAN):");
  disp(mean(M(:)));
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
function [T, ind, d] = sample_image_wavelets(im, J = 10, max_N = 1000)
  ind = get_max_n_local_gradient_extrema_value_indexes(im, max_N);
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

function C = cut_image(im, l, b, offset)
  C = im(offset:offset+b-1, offset:offset+l-1);
endfunction

%%
% params: C, T, fun
function ww = sample_image_series(C, T, fun)
  for i = 1:size(C, 2)
    f = fun(C{i});
    ww{i} = W(f, T);
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
% params: img, ind, d, min_dist = 0
function sample_mesh_image_with_graph(img, ind, d, min_dist = 0)
  sample_mesh_image_with_indices(img, ind);
  hold on;
  for i = 1:size(ind, 1)
    for j = 1:size(ind, 1)
      if (d{i}(j) > min_dist)
        x = [ind(j, 2) ind(i, 2)];
        y = [ind(j, 1) ind(i, 1)];
        plot(x, y, 'r;;');
      endif
    endfor
  endfor
  hold off;
endfunction
