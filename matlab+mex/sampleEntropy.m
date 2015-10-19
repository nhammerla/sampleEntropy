%
% Sample Entropy (mex-version)
%
%   SampEn = sampleEntropy(INPUT, M, R, TAU)
%
%   Calculate sample entropy according to
%   [1] https://en.wikipedia.org/wiki/Sample_entropy. 
%   
%   Arguments:
%       INPUT       Nx1         Input sequence.
%       M           Int         Window-length (or "dimension"). See [1] for
%                               details.
%       R           Double      Tolerance for "similarity". Used as
%                               threshold on cheb. distance [1]. 
%       TAU         Int         Spacing of valid samples (for subsampling).
%                               A value of 1 corresponds to no subsampling,
%                               2 takes every other value, etc.
%       
%   Nils Hammerla '2015 <n.hammerla@gmail.com>