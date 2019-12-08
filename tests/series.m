## function needed that can be used to generate a graph
## where each node is connected to N others...

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

function ww = wavelets_for_image_series(C, T, ind)
  for i = 1:size(C, 2)
    f = function_from_indexes(C{i}, ind);
    ww{i} = W(f, T);
  endfor
endfunction

function [r_index, M] = max_N_values_index(im, N_limit)
  im = abs(im); # make sure we have only positive values
  regional_max = imregionalmax(im).*im;
  sorted_values = wrev(unique(regional_max));

  for k = 1:size(sorted_values, 1)
    maxima_index_im = ismember(im, sorted_values(k));
    [row_index, col_index] = find(maxima_index_im);
    temp_index = [row_index, col_index];
    if (exist("r_index"))
      new_index_size = size(temp_index, 1) + size(r_index, 1);
    else
      new_index_size = size(temp_index, 1);
    endif
    if (new_index_size > N_limit)
      # select a random sample from these indexes
      n_rows = size(temp_index, 1);

      # check if we have already something in the max index
      if (exist("r_index"))
        n_to_choose = N_limit - size(r_index, 1);
      else
        n_to_choose = N_limit;
      endif
      random_rows_index = randperm(n_rows, n_to_choose);
      r_index = [r_index; temp_index(random_rows_index,:)];
      break;
    else
      if (exist("r_index"))
        r_index = [r_index; temp_index];
      else
        r_index = temp_index;
      endif
    endif
  endfor
endfunction

function image_segmentation_fixed_graph_series(C, l, b, img, series, J)
  [chi, lambda] = generate_fixed_graph(l, b);
endfunction

function image_segmentation_series(C, N_limit, img, series, J)
  # 
endfunction

function local_analysis(C, C0, N_limit, img, series, J, type)
  image_type = "png";
  analysis_name = [img "-S" series "-T" type "-J" num2str(J)];
  figure("visible", "off");

  ind = sort_indexes(max_N_values_index(C0, N_limit));
  dist = image_pixel_euclidean_distance_all(ind, true);
  L = laplace_from_distances(dist);

  disp("Calculating eigen values");
  [chi, lambda] = eig(L, "vector");

  plot(lambda);
  print(["-d" image_type], ["~/lambda-" analysis_name]);

  # we are still using the g and h function from the original paper
  disp("Calculating wavelet T matrix");
  T = wavelets(@g, @h, chi, lambda, J);
  ww = wavelets_for_image_series(C, T, ind);

  imagesc(index_to_graph(ind, size(C{1}, 1), size(C{1}, 2)));
  print(["-d" image_type], ["~/graph-" analysis_name]);

  ignore_h_kernel = false;
  print_wavelet_image_compare_series(ww, analysis_name, image_type, ignore_h_kernel);

  print_wavelet_series(ww, analysis_name, image_type, ignore_h_kernel);
endfunction

function difference_series(C, N_limit, img, series, J)
  diffC = C{2} - C{1};
  local_analysis(C, diffC, N_limit, img, series, J, "diff");
endfunction

function gradient_series(C, N_limit, img, series, J)
  gradC = gradient(C{1});
  local_analysis(C, gradC, N_limit, img, series, J, "grad");
endfunction

function [chi_r, lambda_r] = generate_fixed_graph(l, b)
  # Use persistent variables to keep to not recreate them for every test.
  persistent chi;
  persistent lambda;
  persistent L8P;

  if (isempty(L8P))
    L8P = laplace_eight_point(l, b);
  endif

  if (isempty(chi) || isempty(lambda))
    [chi, lambda] = eig(L8P, "vector");
  endif

  chi_r = chi;
  lambda_r = lambda;
endfunction

function run_all_test_series_for(C, img, series, J)
  N_limit = 1000;
  l = 40; b = 25; # to keep consistent with N_limit = 1000.

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
  for j = 1:size(J, 2)
    run_sun_test_series_for(impath, img, "sun-1", J(j));
    #run_sun_test_series_for(impath, img, "sun-2", J(j));
    #run_sun_test_series_for(impath, img, "sun-3", J(j));
    #run_sun_test_series_for(impath, img, "sun-4", J(j));
    #run_sun_test_series_for(impath, img, "sun-5", J(j));
    #run_sun_test_series_for(impath, img, "sun-6", J(j));
  endfor
endfunction

function run_moon_series(impath, img)
  J = [10, 30];
  for j = 1:size(J, 2)
    run_moon_test_series_for(impath, img, "moon-1", J(j));
    #run_moon_test_series_for(impath, img, "moon-2", J(j));
    #run_moon_test_series_for(impath, img, "moon-3", J(j));
    #run_moon_test_series_for(impath, img, "moon-4", J(j));
    #run_moon_test_series_for(impath, img, "moon-5", J(j));
  endfor
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
