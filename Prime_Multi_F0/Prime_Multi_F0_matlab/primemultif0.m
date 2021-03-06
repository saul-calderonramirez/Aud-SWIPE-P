function S = primemultif0(x,fs, plim, dt)
% S = PRIMEMULTIF0(X,Fs,PC,DT) computes 
% pitch scores for the pitch candidates PC
% every DT seconds, based on the signal X
% with sampling frequency FS. PC must be an
% increasing sequence.
timesEx = zeros(1, 5);


tStart = tic;

dlog2p = 1/48;
if length(plim) == 2 % We assume are the limits of the range
    log2pc = ( log2(plim(1)): dlog2p: log2(plim(2)) )';
    pc = 2 .^ log2pc;
else
    pc = plim(:); % We assume PLIM is a SORTED list of candidates
    log2pc = log2(pc);
end

t = [0:dt:length(x)/fs]'; % time
%OJOOOOO, PASA A VECTOR COLUMNA, DEBE FUNCIONAR
%pc = pc(:);
%log2pc = log2(pc);
S = zeros(length(pc),length(t)); % scores
%PORQUE HACE ESE MIN??
limLogWs = round(log2(8*fs./[pc(1) pc(end)]));
ws = 2.^[limLogWs(1):-1:limLogWs(end)]; % window sizes
pO = 8*fs./ws; % optimal pitches for window sizes

 tElapsed = toc(tStart);
 timesEx(1) = tElapsed;
 tags{1} = 'Pitch candidates arrangement time';


% Determine window sizes used by each pitch candidate:
d = 1 + log2pc - log2(8*fs./ws(1));


disp('NUMBER OF WINDOW SIZES');
length(ws)

for i=1:length(ws)
    dn = max(1,round(4*fs/pO(i))); % hop size
    % Zero-pad signal:
    xzp = [zeros(ws(i)/2,1);x(:);zeros(dn+ws(i)/2,1)];
    % Compute spectrum:
    w = hanning(ws(i)); % Hann window 
    %o = woverlap
    o = max(0,round(ws(i) - dn)); % window overlap
    %One fourier transform per column
    %f, frequencies at which the fft was calculated
    tStart = tic;
    
    
    [X,f,ti] = specgram(xzp,ws(i),fs,w,o);
    
    tElapsed = toc(tStart);    
    timesEx(2) = tElapsed + timesEx(2);    
    tags{2} = 'Specgram calculation time';

    
    % Select candidates that use this window size:
    tStart = tic;
    if length(ws) == 1
        j=(1:length(pc))'; k=[];
    elseif i == length(ws)
        j=find(d-i>-1); k=find(d(j)-i<0);
    elseif i == 1
        j=find(d-i<1); k=find(d(j)-i>0);
    else
        j=find(abs(d-i)<1); k=1:length(j);
    end
    % Compute loudness at required frequency range:
    first = find(f>pc(j(1))/4,1,'first');
    f = f(first:end);
    
    tElapsed = toc(tStart); 
    timesEx(3) = tElapsed + timesEx(3);    
    tags{3} = 'Determination of pitch candidates corresponding to ws time';
    
    
    tStart = tic;
    
    %root of the spectogram
    L = sqrt(abs(X(first:end,:)));
    % Compute scores:
    Si = scoresAllCandidates(f,L,pc(j));
    Si2 = interp1(ti,Si',t,'linear',NaN)';
    lambda = d(j(k)) - i;
    mu = ones(size(j));
    mu(k) = 1 - abs(lambda);
    S(j,:) = S(j,:) + repmat(mu,1,size(Si2,2)).*Si2;
    
    tElapsed = toc(tStart); 
    timesEx(4) = tElapsed + timesEx(4);    
    tags{4} = 'Scores of all candidates time';
    
end
% Enhance pitch strength matrix:

tStart = tic;

S = max(0,S);
primesN = primes(pc(end)/pc(1));
for i = primesN;
    primesPcs = i*pc;
    temp = interp1(pc,S, primesPcs,'linear',0);
    temp2 = max(0, temp);
    S = S - max(0, temp2);
end
S = max(0,S);

    
tElapsed = toc(tStart); 
timesEx(5) = tElapsed + timesEx(5);    
tags{5} = 'Post processing time of the scores';


disp('Prime multi F0 finished');

function S = scoresAllCandidates(f,L,pc)
% Create score matrix:
S = zeros(length(pc),size(L,2));
% Define integration regions:
k = ones(1,length(pc)+1);
for j = 1:length(k)-1
    k(j+1) = k(j) - 1 + find(f(k(j):end)>pc(j)/4,1,'first');
end
k = k(2:end);
% Create loudness normalization matrix:
N = (flipud(cumsum(flipud(L))));
for j = 1:length(pc)
    % Normalize loudness:
    n = N(k(j),:);
    n(n==0) = Inf; % to make zero-loudness equal zero after normalization
    NL = L(k(j):end,:)./repmat(n,size(L,1)-k(j)+1,1);
    % Compute score:
    S(j,:) = scoreOneCandidate(f(k(j):end),NL,pc(j));
end
% function S = pitchStrengthAllCandidates( f, L, pc )
% % Create pitch salience matrix
% S = zeros( length(pc), size(L,2) );
% % Define integration regions
% k = ones( 1, length(pc)+1 );
% for j = 1 : length(k)-1
%     k(j+1) = k(j) + find( f(k(j):end) > pc(j)/4, 1, 'first' ) - 1;
% end
% k = k(2:end);
% % Loudness normalization factor
% N = sqrt( flipud( cumsum( flipud(L.*L) ) ) );
% for j = 1 : length(pc)
%     % Normalize loudness
%     n = N(k(j),:);
%     n(n==0) = Inf; % to make zero-loudness equal zero after normalization
%     NL = L(k(j):end,:) ./ repmat( n, size(L,1)-k(j)+1, 1);
%     % Compute pitch strength
%     S(j,:) = pitchStrengthOneCandidate( f(k(j):end), NL, pc(j) );
% end
%  


function S = scoreOneCandidate(f,NL,pc)
n = fix(f(end)/pc-0.75); % number of harmonics
if n==0, S=NaN; return, end
k = zeros(size(f)); % kernel
q = f/pc; % normalized frequency w.r.t. candidate
% Create kernel:
primesN = primes(n);
for i = [1 primes(n)]
    a = abs(q-i);
    % Peak's weigth:
    p = a<0.25; 
    k(p) = cos(2*pi*q(p));
    % Valleys' weights
    v = 0.25<a & a<0.75;
    k(v) = k(v) + cos(2*pi*q(v))/2;
end
normK =  norm( k(k>0) );

k = k / normK; 
S = k'*NL;
pss = 0;

%  
% function S = pitchStrengthOneCandidate( f, NL, pc )
% n = fix( f(end)/pc - 0.75 ); % Number of harmonics
% if n==0, S=0; return, end
% k = zeros( size(f) ); % Kernel
% q = f / pc; % Normalize frequency w.r.t. candidate
% for i = [ 1 primes(n) ] % This is A-SWIPE'. 
%                         % Replace 'primes' with a semicolon for A-SWIPE
%     a = abs( q - i );
%     % Peak's weigth
%     p = a < .25; 
%     k(p) = cos( 2*pi * q(p) );
%     % Valleys' weights
%     v = .25 < a & a < .75;
%     k(v) = k(v) + cos( 2*pi * q(v) ) / 2;
% end
% % Apply envelope
% k = k .* sqrt( 1./f  ); 
% % K+-normalize kernel
% k = k / norm( k(k>0) ); 
% % Compute pitch strength
% S = k' * NL; 
