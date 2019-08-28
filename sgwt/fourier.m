%%
% Graph fourier transform
%
% params: f, chi
function fh = gft(f, chi)
  fh = chi' * f;
endfunction

%%
% Inverse graph fourier transform
%
% params: fh, chi
function f = igft(fh, chi)
  f = chi * fh;
endfunction
