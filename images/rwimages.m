%%
% params: D, name
%
% Write all images in D as name-id.png
function write_all_images(D, name)
  for k = 1:length(D)
    imwrite(D{k}, sprintf('%s-%04d.png', name, k));
  endfor
endfunction

%%
% params: impath
%
% read all *.png images found in the given path
function C = read_all_images(impath)
  ims = dir(fullfile(impath, '*.png'));
  for k = 1:length(ims)
    C{k} = double(rgb2gray(imread(fullfile(impath, ims(k).name))));
  endfor
endfunction

function C = read_all_images_raw(impath)
  ims = dir(fullfile(impath, '*.png'));
  for k = 1:length(ims)
    C{k} = double(imread(fullfile(impath, ims(k).name)));
  endfor
endfunction

