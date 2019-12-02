%%
% params: norm = false
%
% creates a laplace matrix representing a simple icosahedron
function [L, D, A] = laplace_icosahedron(norm = false)
  N = 12;
  D = 5 * eye(N);

  A = [...
  0,1,1,1,1,1,0,0,0,0,0,0;...
  1,0,1,0,0,1,1,1,0,0,0,0;...
  1,1,0,1,0,0,0,1,0,0,0,1;...
  1,0,1,0,1,0,0,0,0,0,1,1;...
  1,0,0,1,0,1,0,0,0,1,1,0;...
  1,1,0,0,1,0,1,0,0,1,0,0;...
  0,1,0,0,0,1,0,1,1,1,0,0;...
  0,1,1,0,0,0,1,0,1,0,0,1;...
  0,0,0,0,0,0,1,1,0,1,1,1;...
  0,0,0,0,1,1,1,0,1,0,1,0;...
  0,0,0,1,1,0,0,0,1,1,0,1;...
  0,0,1,1,0,0,0,1,1,0,1,0;...
  ];

  L = laplace_matrix(D, A, norm);
endfunction

function [L, D, A] = laplace_fully_connected(l, b, norm = false)
  N = l * b;
  A = ones(N, N) - eye(N);
  D = degree_matrix_from_adjacency(A);
  L = laplace_matrix(D, A, norm);
endfunction

function [L, D, A] = laplace_from_distances_fully_connected(l, b, norm = false)
  [r, c] = find(ones(l, b));
  ind = [r c];
  d = image_pixel_euclidean_distance_all(ind);
  [L, D, A] = laplace_from_distances(d, norm);
endfunction

%%
% params: N, norm = false
% same as: laplace_1D(N, true, norm)
%
% creates a laplace matrix representing a circle
function [L, D, A] = laplace_circle(N, norm = false)
  [L, D, A] = laplace_1D(N, true, norm);
endfunction

%%
% params: N, norm = false
% same as: laplace_1D(N, false, norm)
%
% creates a laplace matrix representing a line
function [L, D, A] = laplace_line(N, norm = false)
  [L, D, A] = laplace_1D(N, false, norm);
endfunction

%%
% params: N, connected, norm = false
% 
% N: number of vertices
% connected: are the end vertices connected
% norm: should the laplace matrix be normalized
%
% creates a laplace matrix for a one dimensional graph
% which can optionally be connected to form a circle
function [L, D, A] = laplace_1D(N, connected, norm = false)
  D = 2 * eye(N);
  if (!connected)
    D(1, 1) = 1;
    D(N, N) = 1;
  endif
  A = symetric_band_matrix(N, 1, 1);
  if (connected)
    A(1, N) = 1;
    A(N, 1) = 1;
  endif
  L = laplace_matrix(D, A, norm);
endfunction

%%
% params: N, value, offset
%
% N = matrix size
% value the value to set
% offset offset from the diagonal
%
% create a symetric band matrix containing the value
% with a specified offset
function A = symetric_band_matrix(N, value, offset)
  values = value * ones(N-offset, 1);
  A = diag(values, offset) + diag(values, -offset);
endfunction

%%
% params: l, b, norm = false
%
% l = | longitude, laengengrad |
% b = | latitude, breitengrad |
% create a laplace matrix representing a sphere.
function [L, D, A] = laplace_sphere_grid(l, b, norm = false)
  N = l * b + 2;

  if (l == 1)
    D = 2 * eye(N);
  elseif (l == 2)
    D = 3 * eye(N);
  else
    D = 4 * eye(N);
  endif
  D(1, 1) = l;
  D(N, N) = l;
  A = symetric_band_matrix(N, 1, 1);

  if (l > 1)
    A += symetric_band_matrix(N, l, 1);

    for i = 2:l
      A(1, i) = 1;
      A(i, 1) = 1;
      A(N, N-i+1) = 1;
      A(N-i+1, N) = 1;
    endfor

    for row = 2:l:N-1
      column = row+l-1;
      A(row, column) = 1;
      A(column, row) = 1;
      if (column < N-1)
        A(column, column + 1) = 0;
        A(column + 1, column) = 0;
      endif
    endfor
  endif

  L = laplace_matrix(D, A, norm);
endfunction

%%
% params: A, a, b, value
function A = symetric_set(A, a, b, value)
  A(a, b) = value;
  A(b, a) = value;
endfunction

%%
% params: l, b, norm = false
%
% l = | longitude, laengengrad |
% b = | latitude, breitengrad |
% create a laplace matrix representing a sphere.
function [L, D, A] = laplace_eight_point(l, b, norm = false)
  N = l*b;

  if (l == 1 || b == 1)
    [L, D, A] = laplace_line(N, norm);
    return;
  endif

  % Degree matrix inner points
  D = 8 * eye(N);

  % Degree matrix corner points
  D(1, 1) = 3;
  D(l, l) = 3;
  D(l*(b-1)+1, l*(b-1)+1) = 3;
  D(l*b, l*b) = 3;

  % Degree matrix border points
  % Top and bottom border
  for k = 2:l-1
    D(k, k) = 5;
    D(l*(b-1) + k, l*(b-1) + k) = 5;
  endfor

  % Left and right border
  for k = 1:b-2
    D(k*l+1, k*l+1) = 5;
    D((k+1)*l, (k+1)*l) = 5;
  endfor

  % connect right neighbor
  A = symetric_band_matrix(N, 1, 1);
  % connect bottom neighbor
  A += symetric_band_matrix(N, 1, l);

  if (l - 1 > 1)
    % connect bottom-left neighbor
    A += symetric_band_matrix(N, 1, l - 1);
  endif
  % connect bottom-right neighbor
  A += symetric_band_matrix(N, 1, l + 1);

  if (l-1 > 1)
    for k = 1:b
      right = k*l;

      % remove excessive bottom connections
      A(right, (k-1)*l+1) = 0;
      A((k-1)*l+1, right) = 0;
    endfor  

    for k = 1:b-1
      right = k*l;
      % remove excessive right connections
      A(right, right+1) = 0;
      A(right+1, right) = 0;
    endfor
  endif

  for k = 1:b-2
    right = k*l;
    % remove excessive bottom-right connections
    A(right, right+l+1) = 0;
    A(right+l+1, right) = 0;
  endfor

  L = laplace_matrix(D, A, norm);
endfunction

%%
% params: D, A, norm = false
%
% Generate the laplace matrix from a given degree matrix and
% a given adjacency matrix. Normalize if needed.
function L = laplace_matrix(D, A, norm = false)
  if (norm)
    sqD = D^(-1/2);
    L = eye(size(D, 1)) - sqD * A * sqD;
  else
    L = D - A;
  endif
endfunction

%%
% params: A
%
% A is the adjacency matrix for which the degree matrix
% will be generated
function D = degree_matrix_from_adjacency(A)
  D = diag(sum(A));
endfunction

%%
% params: distances, norm = false
function [L, D, A] = laplace_from_distances(distances, norm = false)
  A(1,1) = 0;
  for i = 1:size(distances, 2)
    for j = 1:size(distances{i}, 2)
      A(i, j) = distances{i}(j);
      A(j, i) = distances{i}(j);
    endfor
  endfor
  D = degree_matrix_from_adjacency(A);
  L = laplace_matrix(D, A, norm);
endfunction
