function k = syntheticSite_kernel(width)
%SYNTHETICSITE_KERNEL Speckle smoothing kernel for apres.syntheticSite.
%
% k = apres.syntheticSite_kernel()        the default 5-bin kernel
% k = apres.syntheticSite_kernel(width)   a Gaussian of width bins
%
% The kernel sets how many bins a speckle grain spans, and so how many
% independent looks a coherence window contains.

if nargin < 1 || isempty(width)
  k = [0.2 0.6 1 0.6 0.2]';
else
  h = ceil(3 * width);
  k = exp(-0.5 * ((-h:h)' / width).^2);
end
end
