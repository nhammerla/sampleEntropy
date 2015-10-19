function d = sampleEntropy(seq, wlen, r, shift)
%
% Sample Entropy (matlab-version)
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
%   NOTE: For long sequences or large data-set use the MEX version of this
%         script, which is approximately 5 times faster.
%
%   Nils Hammerla '2015 <n.hammerla@gmail.com>

% Copyright (c) 2015, Nils Hammerla. n.hammerla@gmail.com
% All rights reserved.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if shift > 1,
    seq = downsample(seq,shift);
end

% allocate space for extracted windows;
D = zeros(length(seq)-(wlen+1), wlen+1);

% extract windows with length wlen+1
for pos=1:length(seq)-wlen-1,
    D(pos,:) = seq(pos:pos+wlen);
end

% initialise
A = 0;
B = 0;

% calculate number of windows with pairwise distance of less than r, for
% two cases:
%   1) B = with windows = 1..wlen 
%   2) A = with windows = 1..wlen+1
for i=1:size(D,1),
    % Chebyshev distance is max(abs(d_ik-d_jk))
    % D(i,i) is 0, but we should not count that.
    % Also D(i,j) is symmetrical (d(i,j)=d(j,i)), therefore we just need to
    % look at D(i+1:end). Effectively we only calculate "half" of the
    % distance matrix. Due to symmetry we can ignore the rest.
    DD = bsxfun(@minus, D(i+1:end,:), D(i,:)); % subtract current window from all future windows.
    DD = abs(DD); % DD now cheb. distance
    
    v1 = max(DD(:,1:end-1),[],2); % maximum along 2nd dim (case 1)
    v2 = max(v1, DD(:,end));      % add last column (case 2)
    
    B = B + sum(v1 < r);
    A = A + sum(v2 < r);
end

% A contains half the matches,
% B contains half the matches. For estimating A/B this doesn't matter
% really.
d = -log(A/B);

end