function Arr = readFromFileAndPlot(path, mPlot)
   fidP = fopen(path);
   Arr = fread(fidP, inf, 'double');
   fclose(fidP);
   if(mPlot == 1)
       figure;
       plot(Arr);
   end
end