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

function d = image_pixel_euclidean_distance(ind, i, j, cutoff = 0.001)
  d = 1/sqrt((ind(i, 1) - ind(j, 1))^2 + (ind(i, 2) - ind(j, 2))^2);
  if (d < 0.001)
    d = 0;
  endif
endfunction

function f = function_from_indexes(im, ind)
  f = 0;
  for i = 1:size(ind, 1)
    f(i) = im(ind(i, 1), ind(i, 2));
  endfor
  f = f';
endfunction
