%% Helpers
% plot psi wavelets as mesh
%
% params: g, t, L, figure_id
function mesh_wavelet_psi(g, t, L)
  mesh_matrix(TfL(g, t, L));
endfunction

function mesh_wavelet_phi(h, L)
  mesh_matrix(TfL(h, 1, L));
endfunction

%%
% params: A, figure_id
function mesh_matrix(A)
  plot_matrix(A, @mesh);
endfunction

function scatter_matrix(A)
  plot_matrix(A, @scatter3);
endfunction

function plot_matrix(A, style_fn)
  [M, N] = size(A);
  x = linspace(0, 2*pi, N);
  y = linspace(0, 2*pi, M);
  [xx, yy] = meshgrid(x, y);
  style_fn(xx, yy, A);
endfunction

function surface_polar_matrix(A)
  [M, N] = size(A);
  theta = linspace(0, 2*pi, N);
  rho = linspace(0, 2*pi, M);
  [x, y] = meshgrid(theta, rho);
  [xx, yy] = pol2cart(x, y);
  surf(xx,yy,A);
endfunction

function v = coord(z, theta, l, f)
  for i = 0:l+1
    v(i+1) = abs(z) * f(i * theta);
  endfor
endfunction

function [ix, iy] = innercoord(xi, theta, b, l)
  for i = 2:b+1
    xc = coord(xi(i), theta, l, @(x) cos(x));
    yc = coord(xi(i), theta, l, @(x) sin(x));
    for j = 1:l+1
      ix(i-1, j) = xc(j);
      iy(i-1, j) = yc(j);
    endfor
  endfor
endfunction

function v = innerchi(chi, i, l)
  v = [chi(i:i+l-1)', chi(i)];
endfunction

function ch = chimatrix(chi, l, b)
  k = 2;
  for i = 2:b+1
    m = innerchi(chi, k, l);
    k += l;
    for j = 1:l+1
      ch(i-1, j) = m(j);
    endfor
  endfor
endfunction

function cc = sphere_matrix(A, l, b)
  cc = [repelem(A(1), l+1)'; chimatrix(A, l, b); repelem(A(size(A, 1)), l+1)'];
endfunction

function surf_sphere(chi, l, b, grid = false)
  [xx, yy, zz] = sphere_coord(l, b);
  cc = [repelem(chi(1), l+1)'; chimatrix(chi, l, b); repelem(chi(size(chi, 1)), l+1)'];
  surf(xx, yy, zz, cc);
  if (!grid)
    shading interp;
  endif
  colorbar();
endfunction

function mesh_sphere(A, l, b)
  [xx, yy, zz] = sphere_coord(l, b);
  cc = sphere_matrix(A, l, b);
  mesh(xx, yy, zz, cc);
  colorbar();
endfunction

function scatter_mesh_sphere(chi, l, b, scale=1.01, s = [])
  [xx, yy, zz] = sphere_coord(l, b);
  mesh(xx, yy, zz, zeros(size(zz)));
  hold on;
  cc = sphere_matrix(chi, l, b);
  scatter3(xx(:)*scale, yy(:)*scale, zz(:)*scale, s, cc(:), "filled");
  hold off;
  colorbar();
endfunction

function scatter_sphere(chi, l, b, scale=1.01, s = [])
  [xx, yy, zz] = sphere_coord(l, b);
  cc = sphere_matrix(chi, l, b);
  scatter3(xx(:)*scale, yy(:)*scale, zz(:)*scale, s, cc(:), "filled");
  colorbar();
endfunction

function [xx, yy, zz] = sphere_coord(l, b)
  phi = pi / (b+1);
  theta = (2*pi) / l;

  for i = 0:b+1;
    z(i+1) = sin(pi/2-i*phi);
    xi(i+1) = cos(pi/2-i*phi);
  endfor
  zz = z'*ones(1, l+1);

  [ix, iy] = innercoord(xi, theta, b, l);
  xx = [repelem(0, l+1)'; ix; repelem(0, l+1)'];
  yy = [repelem(0, l+1)'; iy; repelem(0, l+1)'];
endfunction

