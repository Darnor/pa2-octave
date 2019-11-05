%%
% Calculate the different J values (log10 spaced)
%
% params: J, K = 20, lambda_max = 10, x2 = 2
function [t, lambda_min] = tJ(J, K = 20, lambda_max = 10, x2 = 2)
  lambda_min = lambda_max / K;
  t_min = x2 / lambda_max;
  t_max = x2 / lambda_min;
  t = logspace(log10(t_min), log10(t_max), J);
endfunction

%%
% Calculate the wavelets for the Eigenvectors and Eigenvalues of a Graph-Laplacian
%
% params: g, h, chi, lambda, J, K = 20, lambda_max = 10, x2 = 2
function [T, t, fn_result] = wavelets(g, h, chi, lambda, J, K = 20, lambda_max = 10, x2 = 2)
  must_diag = size(lambda, 1) != 1 && size(lambda, 2) != 1;
  if (must_diag)
    lambda = diag(lambda);
  endif
  lambda_max = max(lambda);
  [t, lambda_min] = tJ(J, K, lambda_max, x2);
  for j = 1:J
    [T{j}, fn_result{j}] = Tf(g, t(j), chi, lambda);
  endfor
  gamma = max(cell2mat(fn_result(:)));
  [T{J+1}, fn_result{J+1}] = Tf(@(x) h(x, gamma, lambda_min), 1, chi, lambda);
endfunction

%%
% Calculate the T-Matrix for a given scaling factor t and the Eigenvectors chi
% and Eigenvalues lambda. The function fn is used to as a Kernel.
%
% params: fn, t, chi, lambda
function [T, fn_result] = Tf(fn, t, chi, lambda)
  fn_result = fn(t*lambda);
  T = chi * diag(fn_result) * chi';
endfunction

%%
% Calculate the wavelet coefficients for the function f
%
% params: f, T
function w = W(f, T, singlevec = false)
  if (singlevec)
    TT = cell2mat(T');
    w = TT * f;
  else
    for i = 1:size(T, 2)
      w{i} = T{i}*f;
    endfor
  endif
endfunction

%%
% Reconstruct the function again from the wavelet coefficients and the given T
%
% params: w, T
function f = iW(w, T, singlevec = false)
  pL = psuedo_inverse(T)
  if (singlevec)
    f = pL * w;
  else
    ww = cell2mat(w(:));
    f = pL * ww;
  endif
endfunction

%%
% Create the pseudo inverse of the cell matrix T
% 
% params: T
function pL = psuedo_inverse(T)
  TT = cell2mat(T');
  pL = (TT'*TT)^(-1)*TT';
endfunction

%%%%%%% scaling function

%%
% g(t*lambda) scaling function.
%
% params: x, alpha = 2, beta = 2, x1 = 1, x2 = 2
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

%%
% Helper for the g(lambda) scaling function.
%
% params: x
function r = s(x)
  r = -5 + 11 * x - 6 * x.^2 + x.^3;
endfunction

%%
% h(lambda) scaling function
%
% params: x, gamma, lambda_min
function r = h(x, gamma, lambda_min)
  r = gamma * e.^(-(x/(0.6*lambda_min)).^4);
endfunction
