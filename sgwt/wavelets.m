%%
% params: fh, chi
function f = igft(fh, chi)
  f = chi * fh;
endfunction

%%
% params: f, chi
function fh = gft(f, chi)
  fh = chi' * f;
endfunction

%%
% function functions
function delta = dirac_delta(x)
  for k=1:size(x,2)
    if (x(k) < 1e-16 && x(k) > -1e-16)
      delta(k) = 1;
    else
      delta(k) = 0;
    endif
  endfor
  delta = delta';
endfunction

function r = f(N)
  x = linspace(0, 2*pi, N);
  r = 10 + sin(x) + 2 * cos(10 * x);
  r = r';
endfunction

%%
% params: J, K, lambda_max, x2 = 2
function [t, lambda_min] = tJ(J, K = 20, lambda_max = 10, x2 = 2)
  lambda_max
  lambda_min = lambda_max / K;
  t_min = x2 / lambda_max;
  t_max = x2 / lambda_min;
  t = logspace(log10(t_min), log10(t_max), J);
endfunction

%%
% params: g, L, J, K = 20, lambda_max = 10, x2 = 2
%function [T, fn_result] = Tft(g, L, J, t)
%  for j = 1:J
%    [T{j}, fn_result{j}] = Tf(g, t(j), L);
%  endfor
%endfunction

function [T, t] = waveletsL(g, h, L, J, K = 20, lambda_max = 10, x2 = 2)
%  [t, lambda_min] = tJ(J, K, lambda_max, x2);
%  [T, fn_result] = Tft(g, L, J, t);
  [chi, lambda] = eig(L);
  lambda = diag(lambda);
  [T, t] = wavelets(g, h, chi, lambda, J, K, lambda_max, x2);
%  for j = 1:J
%    [T{j}, fn_result{j}] = Tf(g, t(j), L);
%  endfor
%  gamma = max(cell2mat(fn_result(:)));
%  [T{J+1}, fn_result{J+1}] = Tf(@(x) h(x, gamma, lambda_min), 1, L);
endfunction

function [T, t] = wavelets(g, h, chi, lambda, J, K = 20, lambda_max = 10, x2 = 2)
  must_diag = size(lambda, 1) != 1 && size(lambda, 2) != 1;
  if (must_diag)
    lambda = diag(lambda);
  endif
  [t, lambda_min] = tJ(J, K, lambda_max, x2);
  for j = 1:J
    [T{j}, fn_result{j}] = Tf(g, t(j), chi, lambda);
  endfor
  gamma = max(cell2mat(fn_result(:)));
  [T{J+1}, fn_result{J+1}] = Tf(@(x) h(x, gamma, lambda_min), 1, chi, lambda);
endfunction

function w = Wf(f, T)
  N = size(T, 2);
  for i = 1:N
    w{i} = T{i} * f;
    figure(i + 50);
    plot(w{i});
  endfor
  %% alternative:
  %TT = cell2mat(T');
  %ww = cell2mat(w(:));
  %w = TT * f;
endfunction

function f = iwavelets(w, T)
  J = size(T, 2) - 1;
  TT = cell2mat(T');
  L = (TT'*TT)^(-1)*TT';
  ww = cell2mat(w(:));
  f = L * ww;
  %f = T{J+1} * w{J+1};
  %for j = 1:J
  %  f += T{j} * w{j};
  %endfor
endfunction

%%
% params: fn, t, L
function [T, fn_result] = TfL(fn, t, L)
  [chi, lambda] = eig(L);
  lambda = diag(lambda);
  [T, fn_result] = Tf(fn, t, chi, lambda);
endfunction

%%
% params: fn, t, chi, lambda
function [T, fn_result] = Tf(fn, t, chi, lambda)
  fn_result = fn(t*lambda);
  T = chi * diag(fn_result) * chi';
  #figure(1);
  #hold on;
  #x = linspace(0, max(lambda), size(lambda, 1));
  #plot(x, fn_result);
  #hold off;
endfunction

%%
% params: f, g, t, L
function r = W(f, g, t, L)
  r = Tf(g, t, L) * f;
endfunction

function r = WPsi(f, psi_t)
  for i = 1:size(psi_t, 2)
    r{i} = psi_t{i} * f;
  endfor
endfunction

%%%%%%% scaling function

%%
% Compute forward transform by explicitly computing eigenvectors and 
% eigenvalues of graph laplacian
%
% see also: sgwt_ftsd(f, g, t, L)
function r = s(x)
  r = -5 + 11 * x - 6 * x.^2 + x.^3;
endfunction

function r = g(x, alpha = 2, beta = 2, x1 = 1, x2 = 2)
  N = size(x);
  for i = 1:N
    if (x(i) < x1)
      r(i) = x1^(-alpha) * x(i)^(alpha);
    elseif (x(i) > x2)
      r(i) = x2^(beta) * x(i)^(-beta);
    else
      r(i) = s(x(i));
    endif
  endfor
  r = r';
endfunction

function r = h(x, gamma, lambda_min)
  gamma
  lambda_min
  r = gamma * e.^(-(x/(0.6*lambda_min)).^4);
endfunction

%%
% params: f, h, gamma, lambda_min, L
function r = S(f, h, gamma, lambda_min, L)
  r = Tf(@(x) h(x, gamma, lambda_min), L) * f;
endfunction

