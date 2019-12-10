%% function needed that can be used to generate a graph
%% where each node is connected to N others...


%%
% sort indexes by Y and X positions
function ind_sorted = sort_indexes(ind)
  [sort_y, orig_i] = sort(ind(:, 2));
  for k = 1:size(orig_i, 1)
    new_row = [ind(orig_i(k), 1) sort_y(k)];
    if (exist("ind_y_sorted"))
      ind_y_sorted = [ind_y_sorted; new_row];
    else
      ind_y_sorted = [new_row];
    endif
  endfor

  [sort_x, orig_j] = sort(ind_y_sorted(:, 1));
  for k = 1:size(orig_j, 1)
    new_row = [sort_x(k) ind_y_sorted(orig_j(k), 2)];
    if (exist("ind_sorted"))
      ind_sorted = [ind_sorted; new_row];
    else
      ind_sorted = [new_row];
    endif
  endfor
endfunction

%%
% Create wavelets for each image and return them.
function ww = wavelets_for_image_series(C, T, ind)
  for i = 1:size(C, 2)
    f = function_from_indexes(C{i}, ind);
    ww{i} = W(f, T);
  endfor
endfunction

function [r_index, new_indexes] = select_max_n(r_index, temp_index, N_limit)
  new_index_size = size(temp_index, 1) + size(r_index, 1);
  if (new_index_size > N_limit)
    % select a random sample from these indexes
    n_rows = size(temp_index, 1);
    n_to_choose = N_limit - size(r_index, 1);
    random_rows_index = randperm(n_rows, n_to_choose);
    temp_index = temp_index(random_rows_index,:);
  endif
  new_indexes = temp_index;
  r_index = [r_index; temp_index];
endfunction

%%
% get the N maximum values given a certain image.
% this function will choose at random if there are
% to many maximas to no introduce bios.
function r_index = max_N_values_index(im, N_limit)
  im = abs(im); % make sure we have only positive values
  regional_max = imregionalmax(im).*im;
  sorted_values = wrev(unique(regional_max));
  r_index = [];
  for k = 1:size(sorted_values, 1)
    maxima_index_im = ismember(im, sorted_values(k));
    [row_index, col_index] = find(maxima_index_im);
    temp_index = [row_index, col_index];
    r_index = select_max_n(r_index, temp_index, N_limit);
    if (size(r_index, 1) == N_limit)
		  disp("Reached index limitation");
      break;
    endif
  endfor
endfunction

function [chi_r, lambda_r] = generate_fixed_eight_point_graph(lparam, bparam)
  % Use persistent variables to prevent the recalculating everything for every test.
  persistent chi;
  persistent lambda;
  persistent L8P;
  persistent l;
  persistent b;

  if (isempty(l))
    l = lparam;
  endif

  if (isempty(b))
    b = bparam;
  endif

  needs_update = (b != bparam || l != lparam);

  if (isempty(L8P) || needs_update)
    disp("Generating Laplace matrix")
    L8P = laplace_eight_point(l, b);
    %L8P = laplace_from_distances_fully_connected(l, b);
  endif

  if (isempty(chi) || isempty(lambda) || needs_update)
    disp("Calculating Eigenvectors and Eigenvalues");
    [chi, lambda] = eig(L8P, "vector");
  endif

  chi_r = chi;
  lambda_r = lambda;
endfunction

%%
% Creates a fixed graph of size lxb including the indexes for each pixel
function [chi, lambda, ind] = prepare_fixed_eight_point_graph(l, b, l_offset, b_offset, l_bound, b_bound)
  [chi, lambda] = generate_fixed_eight_point_graph(l, b);
  for i = 1:l
    for j = 1:b
      l_ind = i + l_offset;
      b_ind = j + b_offset;
      if (l_ind <= l_bound && b_ind <= b_bound)
        if (exist("ind"))
          ind = [ind; b_ind l_ind];
        else
          ind = [b_ind l_ind];
        endif
      endif
    endfor
  endfor
endfunction

function ind = max_N_segregation_index(C, N_limit)
  [N, M] = size(C);
  img_index = [];
  ind = [];

  for k = 1:floor(sqrt(N_limit))
    v = 2*k;
    h = 2*k;
    v_size = max(1, floor(N/v));
    h_size = max(1, floor(M/h));
    R = split_image(C, v, h);
    for i = 1:size(R, 1)
      for j = 1:size(R, 2)
        Q(i, j) = var(R{i, j}(:));
      endfor
    endfor
    [row_index, col_index] = find(Q>median(Q));
    temp_index = [row_index, col_index];
    [img_index, new_indexes] = select_max_n(img_index, temp_index, N_limit);
    ind = [ind; get_point_index(new_indexes, v_size, h_size)];
    if (size(ind, 1) == N_limit)
		  disp("Reached index limitation");
      break;
    endif
  endfor
endfunction

function [dist, ind] = min_distance(direction)
  dist = 0;
  ind = [];
  for k = 1:size(direction, 2)
    if (size(direction{k}, 1) > 0 && direction{k}(1) > dist)
      ind = [direction{k}(2) direction{k}(3)];
      dist = direction{k}(1);
    endif
  endfor
endfunction

%%
% This function is ultra slow but no time for optimizations left.
%
function dist = nearest_direction_distances(ind)
  for i = 1:size(ind, 1)
    N = [];
    NW = [];
    W = [];
    SW = [];
    S = [];
    SE = [];
    E = [];
    NE = [];
  
    current_ind = [ind(i, 1) ind(i, 2)];
    for j = 1:size(ind, 1)
      if (i == j)
        continue;
      endif

      if (ind(j, 1) > current_ind(1, 1) && ind(j, 2) > current_ind(1, 2))
        NE{size(NE, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) == current_ind(1, 1) && ind(j, 2) > current_ind(1, 2))
        N{size(N, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) < current_ind(1, 1) && ind(j, 2) > current_ind(1, 2))
        NW{size(NW, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) < current_ind(1, 1) && ind(j, 2) == current_ind(1, 2))
        W{size(W, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) < current_ind(1, 1) && ind(j, 2) < current_ind(1, 2))
        SW{size(SW, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) == current_ind(1, 1) && ind(j, 2) < current_ind(1, 2))
        S{size(S, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) > current_ind(1, 1) && ind(j, 2) < current_ind(1, 2))
        SE{size(SE, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      elseif (ind(j, 1) > current_ind(1, 1) && ind(j, 2) == current_ind(1, 2))
        E{size(E, 2)+1} = [image_pixel_euclidean_distance(current_ind(1, 1),ind(j, 1), current_ind(1, 2), ind(j, 2)) ind(j, 1) ind(j, 2)];
      endif
    endfor

    [distNE, indNE] = min_distance(NE);
    [distN, indN] = min_distance(N);
    [distNW, indNW] = min_distance(NW);
    [distW, indW] = min_distance(W);
    [distSW, indSW] = min_distance(SW);
    [distS, indS] = min_distance(S);
    [distSE, indSE] = min_distance(SE);
    [distE, indE] = min_distance(E);

    for j = 1:size(ind, 1)
      if (i == j)
        d(j) = 0;
        continue;
      endif
      if (size(indNE, 1) > 0 && ind(j, 1) == indNE(1, 1) && ind(j, 2) == indNE(1, 2))
        d(j) = distNE;
      elseif (size(indN, 1) > 0 && ind(j, 1) == indN(1, 1) && ind(j, 2) == indN(1, 2))
        d(j) = distN;
      elseif (size(indNW, 1) > 0 && ind(j, 1) == indNW(1, 1) && ind(j, 2) == indNW(1, 2))
        d(j) = distNW;
      elseif (size(indW, 1) > 0 && ind(j, 1) == indW(1, 1) && ind(j, 2) == indW(1, 2))
        d(j) = distW;
      elseif (size(indSW, 1) > 0 && ind(j, 1) == indSW(1, 1) && ind(j, 2) == indSW(1, 2))
        d(j) = distSW;
      elseif (size(indS, 1) > 0 && ind(j, 1) == indS(1, 1) && ind(j, 2) == indS(1, 2))
        d(j) = distS;
      elseif (size(indSE, 1) > 0 && ind(j, 1) == indSE(1, 1) && ind(j, 2) == indSE(1, 2))
        d(j) = distSE;
      elseif (size(indE, 1) > 0 && ind(j, 1) == indE(1, 1) && ind(j, 2) == indE(1, 2))
        d(j) = distE;
      else
        d(j) = 0;
      endif
    endfor
    dist{i} = d'; 
  endfor
endfunction

function [chi, lambda] = eig_from_index(ind)
  dist = image_pixel_euclidean_distance_all(ind, true);

  % this is function is just way to slow currently
  %dist = nearest_direction_distances(ind);

  disp("Done calculating distances");
  L = laplace_from_distances(dist);
  [chi, lambda] = eig(L, "vector");
endfunction

function [l_offset, b_offset] = find_offset(C, l, b)
  l_offset = 0;
  b_offset = 0;
  variance = 0;
  base_b_offset = floor(mod(size(C, 1), b)/2)+1;
  base_l_offset = floor(mod(size(C, 2), l)/2)+1;

  for i = 1:size(C, 1)/b
    for j = 1:size(C, 2)/l
      l_start = (j-1)*l+base_l_offset;
      b_start = (i-1)*b+base_b_offset;
      sub_var = var(C(b_start:i*b, l_start:j*l)(:));
      if (sub_var > variance)
        variance = sub_var;
        l_offset = l_start - 1;
        b_offset = b_start - 1;
      endif
    endfor
  endfor
endfunction

function run_analysis(chi, lambda, J, C, ind, img, series, type)
  analysis_name = [img "-S" series "-T" type "-J" num2str(J)];
  image_type = "png";

  figure("visible", "off");
  plot(lambda);
  print(["-d" image_type], ["~/ownCloud/pa2data/" series "/lambda-" analysis_name]);

  % we are still using the g and h function from the original paper
  disp("Calculating wavelet T matrix");
  T = wavelets(@g, @h, chi, lambda, J);
  ww = wavelets_for_image_series(C, T, ind);

  imagesc(index_to_graph(ind, size(C{1}, 1), size(C{1}, 2)));
  print(["-d" image_type], ["~/ownCloud/pa2data/" series "/graph-" analysis_name]);

  ignore_h_kernel = true;
  print_wavelet_image_compare_series(ww, ["~/ownCloud/pa2data/" series "/"], analysis_name, image_type, ignore_h_kernel);

  print_wavelet_series(ww, ["~/ownCloud/pa2data/" series "/"], analysis_name, image_type, ignore_h_kernel);
endfunction

function run_analysis_for_all_J(chi, lambda, J, C, ind, img, series, type)
  for j = 1:size(J, 2)
    run_analysis(chi, lambda, J(j), C, ind, img, series, type);
  endfor
endfunction

function local_analysis(C, C0, N_limit, img, series, J, type)
  ind = sort_indexes(max_N_values_index(C0, N_limit));
  [chi, lambda] = eig_from_index(ind);
  run_analysis_for_all_J(chi, lambda, J, C, ind, img, series, type);
endfunction

%%
% Test series functions
function image_segmentation_fixed_graph_series(C, l, b, img, series, J)
  [l_offset, b_offset] = find_offset(C{1}, l, b);
  [chi, lambda, ind] = prepare_fixed_eight_point_graph(l, b, l_offset, b_offset, size(C{1}, 2), size(C{1}, 1));
  run_analysis_for_all_J(chi, lambda, J, C, ind, img, series, "fixed_graph");
endfunction

function image_segmentation_series(C, N_limit, img, series, J)
  ind = sort_indexes(max_N_segregation_index(C{1}, N_limit));
  [chi, lambda] = eig_from_index(ind);
  run_analysis_for_all_J(chi, lambda, J, C, ind, img, series, "segmentation");
endfunction

function difference_series(C, N_limit, img, series, J)
  local_analysis(C, C{2} - C{1}, N_limit, img, series, J, "diff");
endfunction

function gradient_series(C, N_limit, img, series, J)
  local_analysis(C, gradient(C{1}), N_limit, img, series, J, "grad");
endfunction
% Test series functions end
%%

function run_all_test_series_for(C, img, series, J)
  N_limit = 1000;
  l = 40; b = 25; % to keep consistent with N_limit = 1000.

  disp("Starting segmentation with fixed graph test");
  image_segmentation_fixed_graph_series(C, l, b, img, series, J);

  disp("Starting segmentation test");
  image_segmentation_series(C, N_limit, img, series, J);

  disp("Starting difference test");
  difference_series(C, N_limit, img, series, J);

  disp("Starting gradient test");
  gradient_series(C, N_limit, img, series, J);
endfunction

function run_moon_test_series_for(impath, img, series, J)
  disp(["Starting moon tests. Image: " img ", Series: " series ", J: " num2str(J)]);
  run_all_test_series_for(read_all_images([impath "/" img "/" series "/"]), img, series, J);
  disp("Done with this series.");
endfunction

function run_sun_test_series_for(impath, img, series, J)
  disp(["Starting sun tests. Image: " img ", Series: " series ", J: " num2str(J)]);
  run_all_test_series_for(read_all_images_raw([impath "/" img "/" series "/"]), img, series, J);
  disp("Done with this series.");
endfunction

function run_sun_series(impath, img)
  J = [10, 30];
  run_sun_test_series_for(impath, img, "sun-1", J);
  run_sun_test_series_for(impath, img, "sun-2", J);
  run_sun_test_series_for(impath, img, "sun-3", J);
  run_sun_test_series_for(impath, img, "sun-4", J);
  run_sun_test_series_for(impath, img, "sun-5", J);
  run_sun_test_series_for(impath, img, "sun-6", J);
endfunction

function run_moon_series(impath, img)
  J = [10, 30];
  run_moon_test_series_for(impath, img, "moon-1", J);
  run_moon_test_series_for(impath, img, "moon-2", J);
  run_moon_test_series_for(impath, img, "moon-3", J);
  run_moon_test_series_for(impath, img, "moon-4", J);
  run_moon_test_series_for(impath, img, "moon-5", J);
endfunction

function run_all_moon_series()
  run_moon_series("~/data/20180928", "DSC_0396");

  run_moon_series("~/data/20190125", "DSC_0073");
  run_moon_series("~/data/20190125", "DSC_0074");
  run_moon_series("~/data/20190125", "DSC_0075");
  run_moon_series("~/data/20190125", "DSC_0076");
  run_moon_series("~/data/20190125", "DSC_0077");
  run_moon_series("~/data/20190125", "DSC_0078");
  run_moon_series("~/data/20190125", "DSC_0079");
  run_moon_series("~/data/20190125", "DSC_0080");
  run_moon_series("~/data/20190125", "DSC_0081");
  run_moon_series("~/data/20190125", "DSC_0082");

  run_moon_series("~/data/20190219", "DSC_0083");
  run_moon_series("~/data/20190219", "DSC_0084");
  run_moon_series("~/data/20190219", "DSC_0085");
endfunction

function run_all_sun_series()
  run_sun_series("~/data/20190626", "ASICAP_2019-06-26_18_49_10_643");
  run_sun_series("~/data/20190626", "ASICAP_2019-06-26_18_50_55_924");
  run_sun_series("~/data/20190626", "ASICAP_2019-06-26_18_51_30_884");

  run_sun_series("~/data/20190629", "ASICAP_2019-06-29_07_21_37_333");
  run_sun_series("~/data/20190629", "ASICAP_2019-06-29_07_22_58_258");
  run_sun_series("~/data/20190629", "ASICAP_2019-06-29_07_36_10_397");
  run_sun_series("~/data/20190629", "ASICAP_2019-06-29_07_36_25_606");
  run_sun_series("~/data/20190629", "ASICAP_2019-06-29_07_38_05_483");
endfunction

function run_all_series()
  run_all_moon_series();
  run_all_sun_series();
endfunction
