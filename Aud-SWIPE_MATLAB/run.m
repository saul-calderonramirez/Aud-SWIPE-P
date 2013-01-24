function run
mainDirectory = '../../Muestras/testDefinitivo/';
resultsDirectory = 'results/';
directories = dir(fullfile(mainDirectory));
pmin =  50 ;
pmax =  2000 ;
for m = 1 : length(directories)
    directories(m)
    files = dir(fullfile([mainDirectory directories(m).name '/'], '*.wav' ) );
	for i = 1 : length(files)
        name = files(i).name;
        [x,fs] = wavread([mainDirectory directories(m).name '/' name]);
        [p,t,s] = aswipep( x, fs, name, [pmin pmax], 0.001 );
        data = [ t'; p'];
       
        fid = fopen([resultsDirectory  name '_Result.txt' ],'w');
        fprintf( fid, '%f\t%f\n', data );
        fclose(fid);
        disp([name ' processed.']);
	end
end


beep;