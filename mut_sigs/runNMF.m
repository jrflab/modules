function runNMF( inputFile, outputPrefix, nmfDir, minNumSig, maxNumSig )
% run NMF 
addpath(strcat(nmfDir, '/source/'));
addpath(strcat(nmfDir, '/plotting/'));
clc;

mkdir('temp');

minNumSig = str2num(minNumSig);
maxNumSig = str2num(maxNumSig);

%% Open matlabpool
if ( matlabpool('size') == 0 )
        matlabpool open; % opens the default matlabpool, if it is not already opened
end

%% Define parameters
iterationsPerCore = 100;
stability = zeros(maxNumSig, 1);
reconstructionError = zeros(maxNumSig, 1);
allOutputFile = strcat(outputPrefix, '.mat');

for totalSignatures = minNumSig : maxNumSig
    outputFile = strcat(outputPrefix, '_ts', num2str(totalSignatures), '.mat');

    % Decipher the signatures of mutational processes from catalogues of mutations
    [input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
    decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, ...
    [ outputFile ] );
    % Record the stability and average Frobenius reconstruction error
    stability(totalSignatures-minNumSig+1) = mean(processStabAvg);
    reconstructionError(totalSignatures-minNumSig+1) = norm(input.originalGenomes - processes*exposures, 'fro');
end

%% Plotting the stability and average Frobenius reconstruction error
try %% Some versions of MATLAB plotyy has a bug under linux with -nodisplay -nosplash -nodesktop options
    plotSignatureStabilityAndReconstructionToFile(strcat(outputPrefix, '_stab_reconstruction.png'), minNumSig:maxNumSig, stability, reconstructionError, input);
catch ME
        %% Do not do anything - just ignore the plot in order to save the final output daya
end

%% Saving the data
save(allOutputFile);

quit
end
