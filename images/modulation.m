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

%%
% params: im, fn, l = 'hor' | 'vert', m = 'spline'
function [res, f, xf] = stretch_modulate_image(im, fn, l = 'hor', m = 'spline')
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
        res(j, i) = r(j);
      endfor
    endfor
  else
    x = linspace(1, M, M);
    xf = x+f;
    for i = 1:N
      r = interp1(x, im(i, :), xf, m);
      for j = 1:M
        res(i, j) = r(j);
      endfor
    endfor
  endif
endfunction

%%
% params: im, cutoff_point
function res = cutoff_band_filter(im, cutoff_point)
  for i = 1:size(im, 1)
    for j = 1:size(im, 2)
      if (im(i, j) > cutoff_point || im(i, j) < -cutoff_point)
        res(i, j) = im(i, j);
      else
        res(i, j) = 0;
      endif
    endfor
  endfor
endfunction

%%
% params: im, filter
function r = get_local_maxima_indexes(im, filter)
  M = imregionalmax(filter(im));
  k = 1;
  for i = 1:size(M, 1)
    for j = 1:size(M, 2)
      if (M(i, j) == 1)
        r(k, 1) = i;
        r(k, 2) = j;
        k++;
      endif
    endfor
  endfor
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
    p{k} = p{k} * 1 / max(m);
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
