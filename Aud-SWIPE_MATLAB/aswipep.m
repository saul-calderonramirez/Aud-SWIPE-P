function [p,t,s,pc,S] = aswipep(x,fs,filename,plim,dt,dlog2p,dERBs,woverlap,sTHR)
% ASWIPEP Pitch estimation using A-SWIPE'.
%    P = ASWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR) 
%    estimates the pitch of the vector signal X every DT seconds. The
%    sampling frequency of the signal is Fs (in Hertz). The spectrum is
%    computed using a Hann window with an overlap WOVERLAP between 0 and 1.
%    The spectrum is sampled uniformly in the ERB scale with a step size of
%    DERBS ERBs. The pitch is searched within the range [PMIN PMAX] (in
%    Hertz) with samples distributed every DLOG2P units on a base-2
%    logarithmic scale of Hertz. Pitch estimates with a strength lower than
%    STHR are left undefined.
%    
%    [P,T,S] = ASWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR) 
%    returns the times T at which the pitch was estimated and the pitch 
%    strength S of every pitch estimate.
%
%    P = ASWIPEP(X,Fs) estimates the pitch using the default settings PMIN =
%    30 Hz, PMAX = 5000 Hz, DT = 0.001 s, DLOG2P = 1/48 (48 steps per 
%    octave), DERBS = 0.1 ERBs, WOVERLAP = 0.5, and STHR = -Inf.
%
%    P = ASWIPEP(X,Fs,...,[],...) uses the default setting for the 
%    parameter replaced with the placeholder [].
%
%    REMARKS: (1) For better results, make DLOG2P and DERBS as small as 
%    possible and WOVERLAP as large as possible. However, take into account
%    that the computational complexity of the algorithm is inversely 
%    proportional to DLOG2P, DERBS and 1-WOVERLAP, and that the  default 
%    values have been found empirically to produce good results. Consider 
%    also that the computational complexity is directly proportional to the
%    number of octaves in the pitch search range, and therefore , it is 
%    recommendable to restrict the search range to the expected range of
%    pitch, if any. (2) This code implements A-SWIPE', which uses only the
%    first and prime harmonics of the signal. To convert it into A-SWIPE,
%    which uses all the harmonics of the signal, replace the word
%    PRIMES with a colon (it is located almost at the end of the code).
%    However, this may not be recommendable since A-SWIPE' is reported to 
%    produce on average better results than A-SWIPE (Camacho,2011).
%
%    EXAMPLE: Estimate the pitch of the signal X every 10 ms within the
%    range 75-500 Hz using the default resolution (i.e., 48 steps per
%    octave), sampling the spectrum every 1/20th of ERB, using a window 
%    overlap factor of 50%, and discarding samples with pitch strength 
%    lower than 0.2. Plot the pitch trace.
%       [x,Fs] = wavread(filename);
%       [p,t,s] = aswipep(x,Fs,[75 500],0.01,[],1/20,0.5,0.2);
%       plot(1000*t,p)
%       xlabel('Time (ms)')
%       ylabel('Pitch (Hz)')
%
%    REFERENCES: Camacho, A., (2011) "Using an auditory model to enhance  
%    a sawtooth waveform inspired pitch estimator on signals missing 
%    low-order harmonics," submitted to J. Acoust. Soc. Am.
%
%    NOTE: The method ERBFilters(Fs,F) used here produces the coefficients 
%    of a gammatone filterbank with central frequencies F, given a sampling 
%    frequency FS. The method can be built by modifying the method
%    MakeERBFilters in Malcolm Slaney's Auditory Toolbox. The method
%    ERBFilterBank used here is the one included in that toolbox.
if ~ exist( 'filename', 'var' ) || isempty(filename), filename = 'out'; end
if ~ exist( 'plim', 'var' ) || isempty(plim), plim = [30 5000]; end
if ~ exist( 'dt', 'var' ) || isempty(dt), dt = 0.01; end
if ~ exist( 'dlog2p', 'var' ) || isempty(dlog2p), dlog2p = 1/48; end
if ~ exist( 'dERBs', 'var' ) || isempty(dERBs), dERBs = 0.1; end
if ~ exist( 'woverlap', 'var' ) || isempty(woverlap)
    woverlap = 0.5;
elseif woverlap>1 || woverlap<0
    error('Window overlap must be between 0 and 1.')
end
if ~ exist( 'sTHR', 'var' ) || isempty(sTHR), sTHR = -Inf; end
x = x(:);
t = ( 0: dt: size(x,1)/fs )'; % Times
% Flatten the spectral envelope
b = outmidear( round(fs/10), fs );
[y,zf] = filter( b, 1, x );
y = [ y( fix( (length(b)+1)/2 ) : end ); zf ];
x = y( 1: length(x) );
% Separate groups of consecutive harmonics into channels
f = erbs2hz( 1.5 : hz2erbs(fs/2) )'; % Central frequencies of the channels
fcoefs = ERBFilters(fs,f); % 
X = ERBFilterBank( x, fcoefs )';
% Create new harmonics by applying half-wave rectification
X = resample( X, 2, 1);
X = max( 0, X );
X = resample( X, 1, 2);
% Align the channels
r = round( fs * 4 ./ ( 2*pi*1.019*( 24.7+0.108*f ) ) );
l = size(X,1) - max(r) + 1;
Y = zeros(l,size(X,2));
for k = 1 : size(Y,2)
    Y( :, k ) = X( r(k):r(k)+l-1, k );
end
X = Y;
clear Y;

%Imprimir test 1

 %fid = fopen( fullfile(['./test/1/' filename '.bin' ]),'w');
% fwrite(fid,X,'double');
 %fclose(fid);
% Define pitch candidates
if length(plim) == 2 % We assume are the limits of the range
    log2pc = ( log2(plim(1)): dlog2p: log2(plim(2)) )';
    pc = 2 .^ log2pc;
else
    pc = plim(:); % We assume PLIM is a SORTED list of candidates
    log2pc = log2(pc);
end
S = zeros( length(pc), length(t) ); % Pitch strength matrix
% Determine power-of-two window sizes
logWs = round( log2( 8*fs ./ [ max(-inf,plim(1)) min(inf,plim(end)) ] ) ); 
ws = 2.^( logWs(1): -1: logWs(2) ); % power-of-two window sizes
pO = 8 * fs ./ ws; % Optimal pitches for power-of-two window sizes
% Determine window sizes used by each pitch candidate
d = 1 + log2pc - log2( 8*fs./ws(1) );
% Create ERB-scale uniformly-spaced frequencies (in Hertz)
fERBs = erbs2hz( ( 0: dERBs: hz2erbs(fs/2) )' );
for i = 1 : length(ws)
    % Determine pitch candidates that use this window size
    if length(ws) == 1
        j=( 1:length(pc) )'; k=[];
    elseif i == length(ws)
        j=find(d-i>-1); k=find(d(j)-i<0);
    elseif i==1
        j=find(d-i<1); k=find(d(j)-i>0);
    else
        j=find(abs(d-i)<1); k=1:length(j);
    end
    % Zero pad signal
    dn = max( 1, round( 8*(1-woverlap) * fs / pO(i) ) ); % Hop size
    Xz = [ zeros(ws(i)/2,size(X,2)); X; zeros(dn+ws(i)/2,size(X,2)) ];    
    % Compute specific loudness
    w = repmat( hanning( ws(i) ), 1, size(Xz,2) ); % Hann window
    n = 1 : dn : size(Xz,1)-ws(i)+1; % Centers of the windows
    L = zeros( length(fERBs), length(n) ); % Specific-loudness matrix
    df = fs / ws(i);
    fi = ( 0 : ws(i)-1 )' * df;
    l = ones( 1, length(f)+1 );
    for m = 1 : length(l)-1
        % Compute characteristic frequencies indices
        l(m+1) = l(m) + find( fi(l(m):end) < 1.25*f(m), 1, 'last' ) - 1;
    end
    l = l(2:end);
    W = zeros( length(fi), length(f) );
    for q = 1 : length(f)
        % Compute raised-cosine weights
        W(1:l(q),q) = (1-cos(pi*hz2erbs(fi(1:l(q)))/hz2erbs(f(q))))/2;
    end
    for p = 1 : length(n)
        % Compute specific loudness
        M = abs( fft( w .* Xz( n(p):n(p)+ws(i)-1, : ) ) ); % Magnitudes
        sl = sum( sqrt(W.*M), 2 ); % Specific loudness
        L(:,p) = max( 0, interp1(fi,sl,fERBs,'linear',0) ); % in ERB scale
    end
    % Compute pitch strength
    Si = pitchStrengthAllCandidates( fERBs, L, pc(j) );
    % Interpolate pitch strength at desired times
    if size(Si,2) > 1
        warning off MATLAB:interp1:NaNinY
        Si = interp1( (n-1)/fs, Si', t, 'linear', 0 )';
        warning on MATLAB:interp1:NaNinY
    else
        Si = zeros( length(Si), length(t) );
    end
    % Compute contribution of this window size to pitch strength
    lambda = d( j(k) ) - i;
    mu = ones( size(j) );
    mu(k) = 1 - abs( lambda );
    S(j,:) = S(j,:) + repmat(mu,1,size(Si,2)) .* Si;
end

if length(plim) ~= 2
    % Use only given candidates
    [ s, i ] = max(S,[],1);
    p = pc(i);
    j = find( s < sTHR );
    s(j) = NaN;
    p(j) = NaN;
else
    % Fine-tune the pitch using parabolic interpolation
    p = NaN( size(S,2), 1 );
    s = zeros( size(S,2), 1 );
    for j = 1 : size(S,2)
        [ s(j), i ] = max( S(:,j) );
        if s(j) < sTHR, continue, end
        if i==1 || i==length(pc), p(j)=pc(i); else
            I = i-1 : i+1;
            tc = 1 ./ pc(I);
            ntc = ( tc/tc(2) - 1 ) * 2*pi;
            c = polyfit( ntc, S(I,j), 2 );
            ftc = 1 ./ 2.^( log2(pc(I(1))): 1/12/100: log2(pc(I(3))) );
            nftc = ( ftc/tc(2) - 1 ) * 2*pi;
            [s(j) k] = max( polyval( c, nftc ) );
            p(j) = 2 ^ ( log2(pc(I(1))) + (k-1)/12/100 );
        end
    end
end
 
function S = pitchStrengthAllCandidates( f, L, pc )
% Create pitch salience matrix
S = zeros( length(pc), size(L,2) );
% Define integration regions
k = ones( 1, length(pc)+1 );
for j = 1 : length(k)-1
    k(j+1) = k(j) + find( f(k(j):end) > pc(j)/4, 1, 'first' ) - 1;
end
k = k(2:end);
% Loudness normalization factor
N = sqrt( flipud( cumsum( flipud(L.*L) ) ) );
for j = 1 : length(pc)
    % Normalize loudness
    n = N(k(j),:);
    n(n==0) = Inf; % to make zero-loudness equal zero after normalization
    NL = L(k(j):end,:) ./ repmat( n, size(L,1)-k(j)+1, 1);
    % Compute pitch strength
    S(j,:) = pitchStrengthOneCandidate( f(k(j):end), NL, pc(j) );
end
 
function S = pitchStrengthOneCandidate( f, NL, pc )
n = fix( f(end)/pc - 0.75 ); % Number of harmonics
if n==0, S=0; return, end
k = zeros( size(f) ); % Kernel
q = f / pc; % Normalize frequency w.r.t. candidate
for i = [ 1 primes(n) ] % This is A-SWIPE'. 
                        % Replace 'primes' with a semicolon for A-SWIPE
    a = abs( q - i );
    % Peak's weigth
    p = a < .25; 
    k(p) = cos( 2*pi * q(p) );
    % Valleys' weights
    v = .25 < a & a < .75;
    k(v) = k(v) + cos( 2*pi * q(v) ) / 2;
end
% Apply envelope
k = k .* sqrt( 1./f  ); 
% K+-normalize kernel
k = k / norm( k(k>0) ); 
% Compute pitch strength
S = k' * NL; 
 
function erbs = hz2erbs(hz)
erbs = 21.4 * log10( 1 + hz/229 );

function hz = erbs2hz(erbs)
hz = ( 10 .^ (erbs./21.4) - 1 ) * 229;

