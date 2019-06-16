function plotNMF( prefix, nmfDir, minNumSig, maxNumSig )
% run NMF
addpath(strcat(nmfDir, '/source/'));
addpath(strcat(nmfDir, '/plotting/'));
mkdir('temp');

minNumSig = str2num(minNumSig);
maxNumSig = str2num(maxNumSig);

for totalSignatures = minNumSig : maxNumSig
    tsPrefix = strcat(prefix, '_ts', num2str(totalSignatures));
    inputFile = strcat(tsPrefix, '.mat');
    S = load(inputFile);
    plotSignaturesToFile(tsPrefix, S.processes, S.input, S.allProcesses, S.idx, S.processStabAvg);
    plotSignaturesExposureInSamplesToFile(tsPrefix, S.exposures, S.input);
end

quit
end

