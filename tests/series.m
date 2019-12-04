function image_segmentation_fixed_graph_series(C, l, b)
  [chi, lambda] = generate_fixed_graph(l, b);
endfunction

function image_segmentation_series(C, N_limit)
  # 
endfunction

function difference_series(C, N_limit)
  
endfunction

function gradient_series(C, N_limit)
  
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

function run_all_test_series_for(C)
  N_limit = 1000;
  l = 40; b = 25; # to keep consistent with N_limit = 1000.

  disp("Starting segmentation with fixed graph test");
  image_segmentation_fixed_graph_series(C, l, b);

  disp("Starting segmentation test");
  image_segmentation_series(C, N_limit);

  disp("Starting difference test");
  difference_series(C, N_limit);

  disp("Starting gradient test");
  gradient_series(C, N_limit);
endfunction

function run_moon_test_series_for(impath)
  disp("Starting moon tests");
  run_all_test_series_for(read_all_images(impath));
  disp("Moon series done");
endfunction

function run_sun_test_series_for(impath)
  disp("Starting sun tests");
  run_all_test_series_for(read_all_images_raw(impath));
  disp("Sun series done");
endfunction

function run_sun_series(impath)
  run_sun_test_series_for([impath "/sun-1/"]);
  run_sun_test_series_for([impath "/sun-2/"]);
  run_sun_test_series_for([impath "/sun-3/"]);
  run_sun_test_series_for([impath "/sun-4/"]);
  run_sun_test_series_for([impath "/sun-5/"]);
  run_sun_test_series_for([impath "/sun-6/"]);
endfunction

function run_moon_series(impath)
  run_moon_test_series_for([impath "/moon-1/"]);
  #run_moon_test_series_for([impath "/moon-2/"]);
  #run_moon_test_series_for([impath "/moon-3/"]);
  #run_moon_test_series_for([impath "/moon-4/"]);
  #run_moon_test_series_for([impath "/moon-5/"]);
endfunction

function run_all_moon_series()
  run_moon_series("~/data/20180928/DSC_0396");

  run_moon_series("~/data/20190125/DSC_0073");
  run_moon_series("~/data/20190125/DSC_0074");
  run_moon_series("~/data/20190125/DSC_0075");
  run_moon_series("~/data/20190125/DSC_0076");
  run_moon_series("~/data/20190125/DSC_0077");
  run_moon_series("~/data/20190125/DSC_0078");
  run_moon_series("~/data/20190125/DSC_0079");
  run_moon_series("~/data/20190125/DSC_0080");
  run_moon_series("~/data/20190125/DSC_0081");
  run_moon_series("~/data/20190125/DSC_0082");

  run_moon_series("~/data/20190219/DSC_0083");
  run_moon_series("~/data/20190219/DSC_0084");
  run_moon_series("~/data/20190219/DSC_0085");
endfunction

function run_all_sun_series()
  run_sun_series("~/data/20190626/ASICAP_2019-06-26_18_49_10_643");
  run_sun_series("~/data/20190626/ASICAP_2019-06-26_18_50_55_924");
  run_sun_series("~/data/20190626/ASICAP_2019-06-26_18_51_30_884");

  run_sun_series("~/data/20190629/ASICAP_2019-06-29_07_21_37_333");
  run_sun_series("~/data/20190629/ASICAP_2019-06-29_07_22_58_258");
  run_sun_series("~/data/20190629/ASICAP_2019-06-29_07_36_10_397");
  run_sun_series("~/data/20190629/ASICAP_2019-06-29_07_36_25_606");
  run_sun_series("~/data/20190629/ASICAP_2019-06-29_07_38_05_483");
endfunction

function run_all_series()
  run_all_moon_series();
  run_all_sun_series();
endfunction
