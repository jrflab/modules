function createNMFinput( mutationFile, sampleNameFile, typesFile, cancerType, inputFile)
%create WTSI input 
%   convert mutsig mutation matrix file and sample name file into input for
%   WTSI mutation signature package

originalGenomes = importdata(mutationFile)';

fid = fopen(sampleNameFile);
sampleNames = textscan(fid, '%s');
fclose(fid);
sampleNames = sampleNames{1};

load(typesFile);

save(inputFile, 'originalGenomes', 'subtypes', 'types', 'sampleNames', 'cancerType');
quit
end
