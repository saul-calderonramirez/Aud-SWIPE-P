function mainTest
   [x fs] = wavread('trumpet.wav');
   plim = [150 1175];
   dt = 0.01;
   p0 = readFromFileAndPlot('p0.xlx', 1);
   S = primemultif0(x, fs, plim, dt);
   figure;
   plot(S);
   
end