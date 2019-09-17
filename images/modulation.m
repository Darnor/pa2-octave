pkg load image;

%%
% params: im, fn, l = 'hor' | 'vert'
function [res, f] = multiply_modulate_image(im, fn, l = 'hor')
  [N, M] = size(im);
  isVertical = strcmp(l, 'vert') == 1;
  if (isVertical)
    f = fn(linspace(0,2*pi,N));
  else
    f = fn(linspace(0,2*pi,M));
  endif

  if (isVertical)
    for i = 1:M
      r = im(:, i)'.*f;
      for j = 1:N
        res(j, i) = r(j);
      endfor
    endfor
  else
    for i = 1:N
      r = im(i, :).*f;
      for j = 1:M
        res(i, j) = r(j);
      endfor
    endfor
  endif
endfunction

function S = stretch_modulate_image_series(im, fn, fps, l = 'hor', m = 'spline')
  d = linspace(0, 2*pi, fps);
  for k = 1:fps
    S{k} = stretch_modulate_image(im, @(x) fn(x, d(k)), l, m);
  endfor
endfunction

%%
% params: im, fn, l = 'hor' | 'vert', m = 'spline'
function [S, f, xf] = stretch_modulate_image(im, fn, l = 'hor', m = 'spline')
  [N, M] = size(im);
  isVertical = strcmp(l, 'vert') == 1;
  if (isVertical)
    f = fn(linspace(0,2*pi,N));
  else
    f = fn(linspace(0,2*pi,M));
  endif

  if (isVertical)
    x = linspace(1, N, N);
    xf = x+f;
    for i = 1:M
      r = interp1(x, im(:, i)', xf, m);
      for j = 1:N
        S(j, i) = r(j);
      endfor
    endfor
  else
    x = linspace(1, M, M);
    xf = x+f;
    for i = 1:N
      r = interp1(x, im(i, :), xf, m);
      for j = 1:M
        S(i, j) = r(j);
      endfor
    endfor
  endif
  S(isnan(S)) = 0;
endfunction

%%
% params: im, cutoff_max, cutoff_min
function im = cutoff_band_filter(im, cutoff_max, cutoff_min)
  if (cutoff_max < cutoff_min)
    temp = cutoff_max;
    cutoff_max = cutoff_min;
    cutoff_min = temp;
  endif
  for i = 1:size(im, 1)
    for j = 1:size(im, 2)
      if (im(i, j) > cutoff_max)
        im(i, j) = cutoff_max;
      endif
      if (im(i, j) < cutoff_min)
        im(i, j) = cutoff_min;
      endif
    endfor
  endfor
endfunction

%%
% params: im, cutoff_point
function im = cutoff_min_filter(im, cutoff_point)
  for i = 1:size(im, 1)
    for j = 1:size(im, 2)
      if (im(i, j) < cutoff_point)
        im(i, j) = cutoff_point;
      endif
    endfor
  endfor
endfunction

%%
% params: im, cutoff_point
function im = cutoff_max_filter(im, cutoff_point)
  for i = 1:size(im, 1)
    for j = 1:size(im, 2)
      if (im(i, j) > cutoff_point)
        im(i, j) = cutoff_point;
      endif
    endfor
  endfor
endfunction

%%
% param: im
function r = get_local_extrema_indexes(im)
  reg_min = imregionalmin(im);
  reg_max = imregionalmax(im);
  reg_extremas = reg_min + reg_max;
  k = 1;
  for i = 1:size(reg_extremas, 1)
    for j = 1:size(reg_extremas, 2)
      if (reg_extremas(i, j) >= 1)
        r(k, 1) = i;
        r(k, 2) = j;
        k++;
      endif
    endfor
  endfor
endfunction

%%
% param: im
function r = get_max_n_local_gradient_extrema_value_indexes(im, N)
%  sortedValues = unique(im);               % Unique sorted values
%  maxValues = sortedValues(end-(N-1):end); % Get the 5 largest values
%  maxIndex = ismember(im, maxValues);      % Get a logical index of all values
%                                           % equal to the 5 largest values
  gradient_im = gradient(im);
  reg_min = imregionalmin(gradient_im);
  reg_max = imregionalmax(gradient_im);
  reg_extremas = (reg_min + reg_max).*gradient_im;
  sorted_values = unique(reg_extremas);
  start_idx = 1;
  end_idx = size(sorted_values, 1);
  offset = 0;
  do
    max_values = sorted_values(1+offset:end);
    min_values = sorted_values(1:end-offset);
    maxima_index = ismember(gradient_im, [max_values; min_values]);
    [rows, cols] = find(maxima_index);
    r = [rows cols];
    % cases:
    offset = floor((start_idx + end_idx)/2);
    if (size(r, 1) > N)
      % size(r, 1) > N -> we need less values -> inc start_idx, dec end_idx
      start_idx = offset;
    elseif (size(r, 1) < N)
      % size(r, 1) < N -> we need more values -> dec start_idx, inc end_idx
      end_idx = offset;
    elseif (size(r, 1) == N)
      % size(r, 1) == N -> done
      break;
    endif
  until start_idx >= end_idx;
endfunction

%%
% params: ind, norm = false, filter = @(d) false
function p = image_pixel_euclidean_distance_all(ind, norm = false, filter = @(d) false)
  for i = 1:size(ind, 1)
    for j = 1:size(ind, 1)
      d(j) = image_pixel_euclidean_distance(ind(i, 1), ind(j, 1), ind(i, 2), ind(j, 2));
    endfor
    if (norm)
      m(i) = max(d);
    endif
    for k = 1:size(d, 2)
      if (filter(d(k)))
        d(k) = 0;
      endif
    endfor
    p{i} = d';
  endfor
  for k = 1:size(p)
    p{k} = p{k} / max(m);
  endfor
endfunction

%%
% params: x1, x2, y1, y2
function d = image_pixel_euclidean_distance(x1, x2, y1, y2)
  if (x1 == x2 && y1 == y2)
    % catch diff 0
    d = 0;
  else
    d = 1/sqrt((x1 - x2)^2 + (y1 - y2)^2);
  endif
endfunction

%%
% params: ind, l, b
function G = index_to_graph(ind, l, b)
  G = zeros(l, b);
  for i = 1:size(ind, 1)
    G(ind(i, 1), ind(i, 2)) = 1;
  endfor
endfunction

%%
% params: im, ind
function f = function_from_indexes(im, ind)
  f = 0;
  for i = 1:size(ind, 1)
    f(i) = im(ind(i, 1), ind(i, 2));
  endfor
  f = f';
endfunction
