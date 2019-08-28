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

  if (norm)
    sqD = D^(-1/2);
    L = eye(N) - sqD * A * sqD;
  else
    L = D - A;
  endif
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
  A = symetric_diagonal(N, 1, 1);
  if (connected)
    A(1, N) = 1;
    A(N, 1) = 1;
  endif
  if (norm)
    sqD = D^(-1/2);
    L = eye(N) - sqD * A * sqD;
  else
    L = D - A;
  endif
endfunction

%%
% params: N, value, offset
%
% N = matrix size
% value the value to set
% offset offset from the diagonal
%
% create a symetric diagonal band matrix containing the value
% with a specified offset
function A = symetric_diagonal(N, value, offset)
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
  A = symetric_diagonal(N, 1, 1);

  if (l > 1)
    A += symetric_diagonal(N, l, 1);

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

  if (norm)
    sqD = D^(-1/2);
    L = eye(N) - sqD * A * sqD;
  else
    L = D - A;
  endif
endfunction

