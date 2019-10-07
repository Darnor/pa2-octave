pkg load image;

function fz = focus_measure_brenner(im, n = 2)
  fz = 0;
  [N, M] = size(im);

  for i = 1:N
    for j = 1:M-n
      fz += (im(i, j) - im(i, j+n))^2;
    endfor
  endfor
endfunction

function fz = focus_measure_five_stencil(im, n = 2, m = 2)
  fz = 0;
  [N, M] = size(im);

  for i = 1+n:N-n
    for j = 1+m:M-m
      fz += (4*im(i, j) - im(i, j+m) - im(i, j-m) - im(i+n, j) - im(i-n, j))^2;
    endfor
  endfor
endfunction

function fz = focus_measure_eight_stencil(im, n = 2, m = 2)
  fz = 0;
  [N, M] = size(im);

  for i = 1+n:N-n
    for j = 1+m:M-m
      fz += (...
        8*im(i, j) ...
        - im(i+n, j) ...
        - im(i-n, j) ...
        - im(i, j+m) ...
        - im(i+n, j+m) ...
        - im(i-n, j+m) ...
        - im(i, j-m) ...
        - im(i+n, j-m) ...
        - im(i-n, j-m) ...
        )^2;
    endfor
  endfor
endfunction

function FZ = focus_measure_split(SplitImage, fn)
  FZ = zeros(size(SplitImage));
  for i = 1:size(SplitImage, 1)
    for j = 1:size(SplitImage, 2)
      FZ(i, j) = fn(SplitImage{i, j});
    endfor
  endfor
endfunction

function R = split_image(im, v_splits = 8, h_splits = 8)
  [N, M] = size(im);
  v_size = max(1, floor(N/v_splits));
  h_size = max(1, floor(M/h_splits));

  for j = 1:v_splits
    for i = 1:h_splits
      R{j, i} = im((j-1)*v_size+1:j*v_size, (i-1)*h_size+1:i*h_size);
    endfor
  endfor
endfunction

function ind = get_point_index(r, v, h)
  v_offset = ceil(v / 2);
  h_offset = ceil(h / 2);
  ind = [r(1,1) * v - v_offset, r(1,2) * h - h_offset];
  for i = 2:size(r, 1);
    ind = [ind; r(i,1) * v - v_offset, r(i,2) * h - h_offset];
  endfor
endfunction

function ind = get_index_from_split_image(im, max_N = 1000)
  r = [-1, -1];
  [N, M] = size(im);
  for v = 2:30
    v
    %for h = 8:20
      h = v;
      v_size = max(1, floor(N/v));
      h_size = max(1, floor(M/h));
      R = split_image(im, v, h);
      Q = focus_measure_split(R, @(im) focus_measure_brenner(im));
      [cols, rows] = find(Q>median(Q));
      if (size(cols, 1) + size(r,1) > max_N)
        disp("exceeding max_N");
        disp(size(cols, 1) + size(r,1));
        ind = r;
        return;
      endif
      if (r(1, 1) == r(1, 2) && r(1, 1) == -1)
        r = get_point_index([cols, rows], v_size, h_size);
      else
        r = [r; get_point_index([cols, rows], v_size, h_size)];
      endif
    %endfor
  endfor
  ind = r;
endfunction

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

function R = variance_mean_image_quadrant_internal(im)
  [N, M] = size(im);

  TopLeft = im(1:ceil(N/2), 1:ceil(M/2));
  TopRight = im(1:ceil(N/2), ceil(M/2):M);
  BottomLeft = im(ceil(N/2):N, 1:ceil(M/2));
  BottomRight = im(ceil(N/2):N, ceil(M/2):M);
  
  Var(1, 1) = var(TopLeft(:));
  Var(1, 2) = var(TopRight(:));
  Var(2, 1) = var(BottomLeft(:));
  Var(2, 2) = var(BottomRight(:));

  Std(1, 1) = std(TopLeft(:));
  Std(1, 2) = std(TopRight(:));
  Std(2, 1) = std(BottomLeft(:));
  Std(2, 2) = std(BottomRight(:));

  Mean(1, 1) = mean(TopLeft(:));
  Mean(1, 2) = mean(TopRight(:));
  Mean(2, 1) = mean(BottomLeft(:));
  Mean(2, 2) = mean(BottomRight(:));

  R{1, 1} = Var;
  R{1, 2} = Std;
  R{1, 3} = Mean;
  R{2, 1} = TopLeft;
  R{2, 2} = TopRight;
  R{2, 3} = BottomLeft;
  R{2, 4} = BottomRight;
endfunction

function [ind, R] = variance_mean_image_quadrant(im, max_N = 1024)
  R = variance_mean_image_quadrant_internal(im);
  ind = 0;
  for n = 1:max_N
    m = mean(R{1,2}(:));
    if (R{1,2}(1,1) > m)
      
    endif
    if (R{1,2}(1,2) > m)
      
    endif
    if (R{1,2}(2,1) > m)
      
    endif
    if (R{1,2}(2,2) > m)
      
    endif
  endfor
  ind(cellfun(@isempty, ind)) = [];
%  cutoff = max(starting_variance-(starting_variance/100*percentage), 0);

%  R{k, 1} = Var;
%  R{k, 2} = Mean;

%  if (Var(1, 1) > cutoff)
%    disp("top left");
%    R{k} = variance_mean_image_quadrant(TopLeft, percentage, k+1);
%  endif
%  if (Var(1, 2) > cutoff)
%    disp("top right");
%    R{k} = variance_mean_image_quadrant(TopRight, percentage, k+1);
%  endif
%  if (Var(2, 1) > cutoff)
%    disp("bottom left");
%    R{k} = variance_mean_image_quadrant(BottomLeft, percentage, k+1);
%  endif
%  if (Var(2, 2) > cutoff)
%    disp("bottom right");
%    R{k} = variance_mean_image_quadrant(BottomRight, percentage, k+1);
%  endif
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
    start_idx
    end_idx
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
    else
      m(i) = 1;
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
