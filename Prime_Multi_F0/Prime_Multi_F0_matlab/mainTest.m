function mainTest
   [x fs] = wavread('trumpet.wav');
   plim = [150 3000];
   dt = 0.01;
   %p0 = readFromFileAndPlot('../Prime_Multi_F0_v1/src/PrimeMultiF0/p0.xlx', 1);
   tStart = tic;
   S = primemultif0(x, fs, plim, dt);
   tElapsed = toc(tStart);   
   disp('total execution time')
   tElapsed
   figure;
   plot(S);
   
end