function b = outmidear(n,fs)
% This funtion creates an N-coefficients FIR filter that simulates the
% outer-middle ear. FS is the sampling frequency of the signal to be
% filtered.
% The specification of the filter is taken from a figure in:
% B.R. Glasberg and B.C.J. Moore, A model of loudness applicable to 
% time-varying sounds, J. Audio Eng. Soc. 50: 331-342 (2002)
f = [    0 .02 .05    .1 .2  .5  .6  .7  1  2  3 4.0  5  8   9  10 12  13  14  15] * 1000;
g = [ -inf -39 -19 -12.5 -8  -2  -1   0  0  4  8   8  5 -9 -11 -11 -9 -13 -19 -15]; % gain
m = 10 .^ (g/20);
if 0
    semilogx(f,g,'.-k','linewidth',1);
    set(gca,'xtick',[.02 .05 .1 .2 .5 1 2 5 10]*1000);
    grid on;
    set(gca,'xminorgrid','off');
    set(gca,'xminortick','off');
    xlim([f(1) f(end)]);
end
if (fs/2) > f(end)
    f = [ f (fs/2) ];
    m = [ m 0 ];
    warning(['Nyquist frequency (' num2str(fs/2) ') exceeds filter specifications (0-' num2str(f(end)) ' Hz). Extrapolation towards zero used above ' num2str(f(end)) ' Hz.']);
else
    mend = interp1( f, m, fs/2 );
    i = find( f < fs/2 );
    f = [ f(i) fs/2 ];
    m = [ m(i) mend ];
end
b = fir2( n, f/(fs/2), m, hanning(2*ceil(n/2)+1) );
