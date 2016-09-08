#coding: utf-8
'''

This software is licensed under the following agreement:

JHU Academic Software License Agreement

This license agreement ("License"), effective automatically by and between you (hereinafter referred to as the "Licensee") and The Johns Hopkins University (“JHU”) for use of the software known as CHASM (the “Software”) once and only so long as Licensee complies with the following terms and conditions:

   1. The requirement to acknowledge the copyright of JHU as follows: “Copyright Johns Hopkins University 2010, all rights reserved,” and copyrights of any incorporated third party software as described in the associated documentation.
   2. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
         1. Redistributions of source code must retain the above copyright notice, and these terms and conditions.
         2. Neither the name of JHU nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission from an authorized JHU representative.
         3. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. LICENSEE AGREES TO DEFEND, INDEMNIFY AND HOLD HARMLESS JHU FOR ANY CLAIMS ARISING FROM LICENSEE’S USE OF THE SOFTWARE TO THE FULLEST EXTENT PERMITTED BY LAW.
   3. Only a non-exclusive, nontransferable license is granted to Licensee to use the Software for academic, non-profit, or government-sponsored research purposes. Use of the Software under this License is restricted to non-commercial purposes. Commercial use of the Software requires a separately executed written license agreement.
   4. Licensee agrees that it will use the Software, and any modifications, improvements, or derivatives to the Software that the Licensee may create (collectively, "Improvements") solely for non-commercial purposes and shall not distribute or transfer the Software or Improvements to any person or third parties without requiring that such third parties adhere to the terms of this License.
   5. Licensee acknowledges that JHU holds copyright in the Software or portions of the Software, and that the Software incorporates third party software which may be subject to additional terms and conditions. Licensee agrees that any Improvements made by Licensee shall be subject to the same terms and conditions as the Software.
   6. Licensee agrees that any publication of results obtained with the Software will acknowledge its use by an appropriate citation.
   7. Licensee’s rights under this License terminate automatically without notice from JHU if you fail to comply with any term(s) of this License.
   8. This License shall be governed by the laws of the State of Maryland, excluding the application of its conflicts of law rules. Licensee agrees that any dispute shall be appropriate only in the state and federal courts located within the State of Maryland.

Created May 17, 2010

Contact: chasm-users@lists.johnshopkins.edu
'''

import re, string
import sys
import os.path

class DataSet:
    '''
    Dataset Object to keep track of files and classes
    Gets loaded when the XML is parsed

    MutationFileParser Function
    Initialize with location of the file
    Check that a line is valid using regular expression
    Return the parsed line: Transcript, Position, Wildtype, Mutation

    Verify if the mutation is valid
    MPO = MutationParserObject("some directory path")
    for Transcript, Position, Wildtype, Mutation in MPO.read():
        if dbhelper.validate(Transcript, Position, Wildtype, Mutation):

    If there is only 1 class, then by default no classes will be specified
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.dataset = {} # {Class: [list of files]}

    def addFile(self, datafile, dataclass=0):
        '''
        If the specified file exists
        If class exists, add datafile to existing class
        If class does not exist, creates a new class and adds datafile to that class
        Note: datafile here is the full path
        '''
        if os.path.isfile(datafile):
            if dataclass not in self.dataset.keys():
                self.dataset[dataclass] = datafile
            else:
                raise Exception("Two or more files assigned to same class. Please submit 1 file per class\n")
        else:
            raise IOError("Data file not found")

    def getClasses(self):
        '''
        Function to output a sorted list of classes
        '''
        classes = self.dataset.keys()
        classes.sort()
        return classes

    def getFiles(self, dataclass):
        '''
        Function to retrieve files for a particular class
        '''
        if dataclass in self.dataset.keys():
            return self.dataset.get(dataclass)
        else:
            return None

    def parseTranscriptFile(self, filename):
        '''
        Read missense mutations into a list
        '''
        mutations = []
        for mutation in MutationReader(filename):
            if mutation != None:
                if self.__isValidTranscript(mutation, filename):
                    mutations.append(mutation)
        return mutations

    def __isValidTranscript(self, mut, filename):
        '''
        Check whether mutations are properly formatted
        '''
        status = False
        # Allow for inclusion of uniq ID but do not require it
        m1 = re.search("([A-Za-z0-9\\_\\-\\.]*)\t([A-Z])([0-9]*)([A-Z])",mut.line)
        m2 = re.search("([A-Za-z0-9\\_\\-\\.]*)\t([A-Za-z0-9\\_\\-\\.]*)\t([A-Z])([0-9]*)([A-Z])",mut.line)
        if m1 is not None or m2 is not None:
            status = True
        else:
            sys.stderr.write("Warning: Line " +str(mut.index)+ " is not formatted properly. The line will be skipped. \n")
        return status

    def parseGenomicFile(self, filename):
        '''
        Read genomic single nucleotide variants into a list
        '''
        snvs = []
        for snv in SNVReader(filename):
            if self.__isValidGenomic(snv, filename):
                snvs.append(snv)
        return snvs

    def __isValidGenomic(self, snv, filename):
        '''
        Check whether genomic single nucleotide variants are properly formatted
        '''
        status = False
        # Allow for inclusion of uniq ID but do not require it
        m1 = re.search("([A-Za-z0-9]*)\t([0-9]*)\t([-\+])\t([A-Z])\t([A-Z])$", snv.line)
        m2 = re.search("([A-Za-z0-9]*)\t([A-Za-z0-9]*)\t([0-9]*)\t([-\+])\t([A-Z])\t([A-Z])$", snv.line)
        if m1 is not None or m2 is not None:
            status = True
        else:
            sys.stderr.write("Warning: Line " +str(snv.index)+ " is not formatted properly. The line will be skipped. \n")
        return status

    def getMutations(self): # Currently only used for testing
        '''
        Generator function that parses all files and returns
        mutations and class label
        '''
        for dataclass in self.dataset.keys():
            currentdata = open(self.dataset.get(dataclass),"r")
            line = currentdata.readline()
            linenumber = 1
            while line != "":
                m = re.search("([A-Za-z0-9\\_\\-\\.]*)\t([A-Z])([0-9]*)([A-Z])",line)
                # If this mutation is contained within the database
                if m is not None:
                    yield m.group(1), int(m.group(3)), m.group(2), m.group(4), dataclass
                else:
                    #if len(self.dataset.keys()) == 1:
                    #    sys.stderr.write("Class " + str(dataclass) + " Line " + str(linenumber) + ":\n")
                    #else:
                    #    sys.stderr.write("Line " + str(linenumber) + ":\n")
                    sys.stderr.write("\t\"" + line + "\" is not a properly formatted mutation\n")
                line = currentdata.readline()
                linenumber += 1
            currentdata.close()

class MutationReader(object):
    """
    Read file with transcript and
    amino acid substitutions
    """
    def __init__(self, filename):
        self.fh = file(filename, 'r')

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.fh.readline().strip()
            if line == '':
                self.fh.close()
                raise StopIteration
            line = line.replace(' ', '\t')
            return Mutation(line)

class Mutation(object):
    """
    Stores missense mutation
    """
    index = 1

    def __init__(self, line):
        # Format of line:
        # uid <tab> mRNA accession <tab> codon change <tab> AA position <tab> reference AA <tab> alternate AA <tab> AA len
        #
        self.index = Mutation.index
        Mutation.index += 1
        self.line = line
        self.genomicstatus = 0
        if line.find("\t") != -1:
            toks = self.line.split('\t')
            # If no UID is supplied and if transcript coordinates are used
            if len(toks) == 2:
                self.uid = str(self.index)
                self.acc = toks[0]
                self.ref_aa = toks[1][0]
                self.aa_pos = toks[1][1:-1]
                self.alt_aa = toks[1][-1]
            elif len(toks) == 3:
                self.uid = toks[0]
                self.acc = toks[1]
                self.ref_aa = toks[2][0]
                self.aa_pos = toks[2][1:-1]
                self.alt_aa = toks[2][-1]
            elif len(toks) == 4:
                self.uid = toks[0]
                self.acc = toks[1]
                self.ref_aa = toks[2][0]
                self.aa_pos = toks[2][1:-1]
                self.alt_aa = toks[2][-1]
            #elif len(toks) == 7:
            #    self.uid = toks[0]
            #    self.acc = toks[1]
            #    self.codon_change = ''
            #    self.aa_pos = toks[2]
            #    self.ref_aa = toks[3]
            #    self.alt_aa = toks[4]
            #    self.aa_len = toks[5]
            #    self.codon_strand = toks[6]
            elif len(toks) == 8:
                self.uid = toks[0]
                self.acc = toks[1]
                self.codon_change = toks[2]
                self.aa_pos = toks[3]
                self.ref_aa = toks[4]
                self.alt_aa = toks[5]
                self.aa_len = toks[6]
                self.codon_strand = toks[7]
            self.UID = None
        else:
            raise Exception("Expected a tab delimited file")

    def storeGenomicCoords(self, genomic):
        """
        If generated from a genomic event by snvGetTranscriptList,
        store genomic coordinates
        """
        self.chrom = genomic.chrom
        #self.start = genomic.start
        self.end = genomic.end
        self.strand = genomic.strand
        self.refbase = genomic.ref
        self.altbase = genomic.alt
        self.genomicstatus = 1

    def __repr__(self):
        if hasattr(self, 'codon_change'):
            return "\t".join([str(self.uid), self.acc, self.codon_change, self.aa_pos, self.ref_aa, self.alt_aa, str(self.UID), self.chrom, self.end, self.strand, self.refbase, self.altbase])
        else:
            return "\t".join([str(self.uid), self.acc, self.aa_pos, self.ref_aa, self.alt_aa, str(self.UID)])


class SNVReader(object):
    """
    Read a file with genomic coordinates
    """
    def __init__(self, filename):
        self.fh = file(filename, 'r')

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.fh.readline()
            if line == "":
                self.fh.close()
                raise StopIteration
            line=line.strip().replace(' ','\t')
            return SNV(line)


class SNV(object):
    """
    Store genomic coordinate
    """
    index =  1

    def __init__(self, line):
        self.index = SNV.index
        SNV.index += 1
        self.line = line
        row = self.line.split("\t")
        if len(row) > 0:
            # If user did not assign uid, store
            # uid as row index
            i = 0
            if len(row) == 5:
                self.uid = str(self.index)
                i = 1
            # otherwise, store user assigned
            # uid
            elif len(row) == 6:
                self.uid = row[0]
            self.chrom = row[1-i]
            #self.start = row[2-i]
            self.end = row[2-i]
            self.strand = row[3-i]
            self.ref = row[4-i]
            self.alt = row[5-i]
        else:
            raise Exception("Expected a tab delimited file")

if __name__ == "__main__":
    '''
    Some test code
    '''
    ds = DataSet()
    ds.addFile("../example/Passengers.ids", "passenger")
    ds.addFile("../example/Drivers.ids", "driver")
    print ds.getFiles("passenger")
    print ds.getFiles("driver")
    for transcript, position, wildtype, mutation, dataclass in ds.getMutations():
        print dataclass
        print transcript
        print position
        print wildtype
        print mutation

