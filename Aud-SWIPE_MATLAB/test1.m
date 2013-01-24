function test1
   testsJIAE
   
 function testsJIAE
    directories = {'test/resultsTrumpetC/', 'test/resultsTrumpetC_bonito'};
    fileExt = '*.wav_Result.txt';
    fileExtML = '.*'
    directoryML = directories{1};   
    filesML = dir( fullfile( directoryML ) );
    filesML = filesML(3: size(filesML, 1));
    nFilesML = length(filesML);    
    directoryP = directories{2};   
    filesP = dir( fullfile( directoryP, fileExt ) );
    nFilesP = length(filesP);
    fid = fopen('test/Final_Results.txt','w+');
    for i = 1 : nFilesML        
        nameML = filesML(i).name;
        nameP = filesP(i).name;
        fullnameML = fullfile( directoryML, nameML );
        fullnameP = fullfile( directoryP, nameP );
        P = leerAlturasArchivo(fullnameP);
        %c results have one sample more     
        ML = leerAlturasArchivo(fullnameML); 
        if(length(ML) < length(P))
            P = P(1:size(P,2)-1);   
        end
        if(length(P) == length(ML))
            errores(i) = jIAE(ML, P);
            fprintf( fid, '%s\t%f\t\n', nameML, errores(i) );
        else
            NoHizotest = 1
        end
    end
   
    res = mean(errores)
    desv = std(errores);
    fprintf( fid, '%s\t%f\t\n', 'Promedio de error entre matlab y C (Test 1)', res );
    fprintf( fid, '%s\t%f\t\n', 'Desviacion Estandar de error matlab y C (Test 1)', desv );
fclose(fid);
disp('Done.');
beep;
    
    



function plotFile(fileName)
    file1 = ['test/ResultsTrumpetC/' fileName]
    file2 = ['test/resultsTrumpetML/' fileName]
    values1 = leerAlturasArchivo(file1);
    values2 = leerAlturasArchivo(file2);
    figure;
    plot(values1);
    title([fileName 'fs/10']);
    figure;    
    plot(values2);
    title([fileName 'fs/1000']);    
    
    
function error = jIAE(ML, C)
    error = mean(abs(double(ML) - double(C))./double(ML), 2);
 

function testsAC
    directories = {'test/resultsTrumpetML/', 'test/ResultsTrumpetC/'};
    fileExt = '*.wav_Result.txt';
    fileExtML = '.*'
    directoryML = directories{1};   
    filesML = dir( fullfile( directoryML ) );
    filesML = filesML(3: size(filesML, 1));
    nFilesML = length(filesML);    
    directoryP = directories{2};   
    filesP = dir( fullfile( directoryP, fileExt ) );
    nFilesP = length(filesP);
    fid = fopen('test/Final_Results.txt','w+');


    dev  = 0;   
    ind = 1;

    for i = 1 : nFilesML        
        nameML = filesML(i).name;
        nameP = filesP(i).name;
        fullnameML = fullfile( directoryML, nameML );
        fullnameP = fullfile( directoryP, nameP );
        P = leerAlturasArchivo(fullnameP);
        %c results have one sample more
        P = P(1:size(P,2)-1);
        P = log2(P);
        ML = leerAlturasArchivo(fullnameML); 
        ML = log2(ML);
        ML1 = ML - mean(ML);
        ML = ML1./(norm(ML1));
        P1 = P - mean(P);
        P = P1./norm(P1);    
        sizeML = size(ML, 2);
        sizeP = size(P, 2);
        if(sizeML == sizeP)
            C = xcorr(ML, P, 'coeff');
            Ma = max(C);           
            if(isnan(Ma) == 0)                
                dev(ind) = Ma;                
                ind = ind + 1;
                fprintf( fid, '%s\t%f\t\n', nameML, Ma );

                nameML = nameML(1: size(nameML, 2)-8);
                nameP = nameP(1: size(nameP, 2)- 15);
                 if( nameML == nameP)
                    disp([nameML 'processed']); 
                 else
                    disp([nameML 'Not processed. ' num2str(Ma) ' X' ]); 
                 end
            else
                disp([nameML ' Error in Ma. ' num2str(Ma) ' check' ]);
            end
        end
    
    end
   
    res = mean(dev)
    desv = std(dev);
    fprintf( fid, '%s\t%f\t\n', 'Promedio de correlacion maxima entre matlab y C (Test 1)', res );
    fprintf( fid, '%s\t%f\t\n', 'Desviacion Estandar de correlacion maxima entre matlab y C (Test 1)', desv );
fclose(fid);
disp('Done.');
beep;


function valores = leerAlturasArchivo(nombreArchivo)
    fid = fopen(nombreArchivo);  
    i = 1;
    tline = fgetl(fid);
    while(ischar(tline))    
        [remain token] = strtok(tline);
        valores(i) = str2double(token);
        if(isnan( valores(i))|| valores(i) == 150)
            valores(i) = 1;
        end
        tline = fgetl(fid);
        i = i + 1;
    end    
    fclose(fid);