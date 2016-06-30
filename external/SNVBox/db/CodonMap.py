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

Created Nov 10, 2011

Contact: chasm-users@lists.johnshopkins.edu
'''

import math
import sys
import random

# external imports
from DataSet import Mutation

# vars
rcDict = {"A":"T","T":"A","C":"G","G":"C","N":"N","a":"t","t":"a","g":"c","c":"g","n":"n"}
codonTable = {"ATG":"M", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", "TGT":"C", "TGC":"C", "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "TTT":"F", "TTC":"F", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", "CAT":"H", "CAC":"H", "ATT":"I", "ATC":"I", "ATA":"I", "AAA":"K", "AAG":"K", "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "AAT":"N", "AAC":"N", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R", "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "TGG":"W", "TAT":"Y", "TAC":"Y", "TGA":"*", "TAA":"*", "TAG":"*","AUG":"M", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A", "UGU":"C", "UGC":"C", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "UUU":"F", "UUC":"F", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G", "CAU":"H", "CAC":"H", "AUU":"I", "AUC":"I", "AUA":"I", "AAA":"K", "AAG":"K", "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "AAU":"N", "AAC":"N", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "UGG":"W", "UAU":"Y", "UAC":"Y", "UGA":"*", "UAA":"*", "UAG":"*"}

base_rev_dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'X':'X', '.':'.', 'N':'N', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'x':'x', 'n':'n'}

class CodonConverter(object):
    """
    Convert genomic coordinates and a nucleotide substitution into 
    a transcript and amino acid substitution.

    Note: Only missense variation will be returned unless missenseOnly is set to False.
    """

    def __init__(self, dbconn, missenseOnly=False, pickOneTranscript=False):
        self.conn = dbconn
        self.missenseOnly = missenseOnly
        self.pickOneTranscript = pickOneTranscript
        self.codon_display_len_cutoff = 40
        
    def convert(self, genomic):
        """Get codon, make mutation, return transcript + aa subst if mutation is missense """

        mutations = []
        best_mutation = None

        (codons, best_codon) = self.__get_codons(genomic)

        if best_codon != None:
            best_mutation = Mutation('\t'.join([genomic.uid, \
                                                best_codon.mrnaAcc, \
                                                best_codon.codon_change, \
                                                str(best_codon.codonNo), \
                                                best_codon.ref_aa, \
                                                best_codon.alt_aa, \
                                                str(best_codon.aa_len), \
                                                best_codon.strand]))
            best_mutation.storeGenomicCoords(genomic)
            best_mutation.sequence_ontology = best_codon.sequence_ontology

        for codon in codons:

            mutation = None

            if codon != None:

                # Checks if the reference base matches.
                if codon.sequence_ontology in ['MS', 'SY', 'SG', 'SL'] and self.__matches(genomic, codon) == False:
                    #sys.stderr.write(str(genomic.uid) + "\t" + "reference_base_mismatch\n")
                    continue

                if (self.missenseOnly and codon.sequence_ontology == 'MS') or self.missenseOnly == False:
                    mutation = Mutation("\t".join([genomic.uid, \
                                        codon.mrnaAcc, \
                                        codon.codon_change, \
                                        str(codon.codonNo), \
                                        codon.ref_aa, \
                                        codon.alt_aa, \
                                        str(codon.aa_len), \
                                        codon.strand]))
                    mutation.storeGenomicCoords(genomic)
                    mutation.sequence_ontology = codon.sequence_ontology
            else:
                #sys.stderr.write(str(genomic.uid) + "\t" + "no_codon_mapping\n")
                pass

            if mutation != None:
                mutations.append(mutation)

        #print '### at the end of CodonMap.convert. mutations=', type(mutations)

        return (mutations, best_mutation)

    def __get_reverse_seq (self, seq):
        new_seq = ''
        for i in xrange(len(seq)-1, -1, -1):
            new_seq += base_rev_dic[seq[i]]
        return new_seq

    def __get_ref_codon_bases (self, codon, genomic):
        if codon.strand == genomic.strand:
            start_portion = codon.nucs[:codon.codon_pos].lower()
            middle_portion = genomic.ref.upper()
            codon_display = start_portion + middle_portion
            len_codon_display = len(codon_display)
            len_codon_display_3_div_rem = len_codon_display % 3
            if len_codon_display <= 3:
                last_portion = codon.nucs[len_codon_display:].lower()
            else:
                if len_codon_display_3_div_rem == 0:
                    last_portion = ''
                else:
                    codon_no_increase = len_codon_display / 3
                    last_codon_bases = CodonTblRow.get_codon_bases(self.conn, codon.uid, codon.codonNo + codon_no_increase)
                    if last_codon_bases == None:
                        last_portion = '???'
                    else:
                        last_portion = last_codon_bases[len_codon_display_3_div_rem:].lower()
            codon_display += last_portion
            return codon_display, codon.codonNo, start_portion, middle_portion, last_portion
        else:
            rev_change_bases = self.__get_reverse_seq(genomic.ref)
            middle_portion = rev_change_bases.upper()
            last_portion = codon.nucs[codon.codon_pos + 1:].lower()
            codon_display = middle_portion + last_portion
            len_codon_display = len(codon_display)
            len_codon_display_3_div_rem = len_codon_display % 3
            if len_codon_display <= 3:
                first_codon_no = codon.codonNo
                if len_codon_display_3_div_rem == 0:
                    first_portion = ''
                else:
                    first_portion = codon.nucs[:3 - len_codon_display_3_div_rem].lower()
            else:
                if len_codon_display_3_div_rem == 0:
                    codon_no_decrease = len_codon_display / 3 - 1
                    first_codon_no = codon.codonNo - codon_no_decrease
                    first_portion = ''
                else:
                    codon_no_decrease = len_codon_display / 3
                    first_codon_no = codon.codonNo - codon_no_decrease
                    first_codon_bases = CodonTblRow.get_codon_bases(self.conn, codon.uid, first_codon_no)
                    first_portion = first_codon_bases[:3 - len_codon_display_3_div_rem].lower()
            codon_display = first_portion + codon_display
            return codon_display, first_codon_no, first_portion, middle_portion, last_portion

    def __set_ref_alt_codons(self, codon, genomic):
        codon.ref_codon, codon.codonNo, ref_codon_first_portion, ref_codon_middle_portion, ref_codon_last_portion = self.__get_ref_codon_bases(codon, genomic)

        if codon.strand != genomic.strand:
            alt_seq = self.__get_reverse_seq(genomic.alt)
        else:
            alt_seq = genomic.alt
        codon.alt_codon = ref_codon_first_portion + alt_seq.upper() + ref_codon_last_portion

        if len(codon.ref_codon) > self.codon_display_len_cutoff:
            codon.ref_codon = str(len(codon.ref_codon)) + 'nt'
        if len(codon.alt_codon) > self.codon_display_len_cutoff:
            codon.alt_codon = str(len(codon.alt_codon)) + 'nt'
        codon.codon_change = codon.ref_codon + '>' + codon.alt_codon

    def __set_ref_alt_aas(self, codon, genomic):
        ref_codon_len = len(codon.ref_codon)
        alt_codon_len = len(codon.alt_codon)
        if ref_codon_len == 3 and alt_codon_len == 3:
            codon.ref_aa = self.__translate_codon_bases(codon.ref_codon.upper())
            codon.alt_aa = self.__translate_codon_bases(codon.alt_codon.upper())
        else:
            codon.ref_aa = ''
            codon.alt_aa = ''

    def __set_sequence_ontology(self, codon, genomic):
        ref_seq_len = len(genomic.ref)
        alt_seq_len = len(genomic.alt)
        ref_codon_len = len(codon.ref_codon)
        alt_codon_len = len(codon.alt_codon)
        if ref_seq_len == 1 and alt_seq_len == 1:
            ref_aa = codon.ref_aa
            alt_aa = codon.alt_aa
            if ref_aa == '*':
                if alt_aa == '*':
                    codon.sequence_ontology = 'SY' #'synonymous_variant'
                else:
                    codon.sequence_ontology = 'SL' #'stop_lost'
            elif alt_aa == '*':
                if ref_aa == '*':
                    codon.sequence_ontology = 'SY' #'synonymous_variant'
                else:
                    codon.sequence_ontology = 'SG' #'stop_gained'
            else:
                if ref_aa == alt_aa:
                    codon.sequence_ontology = 'SY' #'synonymous_variant'
                else:
                    codon.sequence_ontology = 'MS' #'missense_variant'
        else:
            if ref_seq_len == 1 and alt_seq_len > 1:
                if (alt_seq_len - ref_seq_len) % 3 == 0:
                    codon.sequence_ontology = 'II' # Inframe insertion
                else:
                    codon.sequence_ontology = 'FI' # Frameshift variant
            elif alt_seq_len == 1 and ref_seq_len > 1:
                if (ref_seq_len - alt_seq_len) % 3 == 0:
                    codon.sequence_ontology = 'ID' # Inframe deletion
                else:
                    codon.sequence_ontology = 'FD' # Frameshift variant
            else:
                codon.sequence_ontology = 'CS'

    def __pick_representative_codon(self, codons):
        if len(codons) == 0:
            return None

        best_transcript_codon = codons[0]
        for codon in codons:
            transcript_mrna_acc_type = codon.mrna_acc_type
            best_transcript_mrna_acc_type = best_transcript_codon.mrna_acc_type
            if transcript_mrna_acc_type == 'RefSeq':
                if best_transcript_mrna_acc_type != 'RefSeq':
                    best_transcript_codon = codon
                else:
                    if codon.aa_len > best_transcript_codon.aa_len:
                        best_transcript_codon = codon
            elif transcript_mrna_acc_type == 'Ensembl':
                if best_transcript_mrna_acc_type == 'RefSeq':
                    pass
                elif best_transcript_mrna_acc_type == 'Ensembl':
                    if codon.aa_len > best_transcript_codon.aa_len:
                        best_transcript_codon = codon
                elif best_transcript_mrna_acc_type == 'CCDS':
                    best_transcript_codon = codon
            elif transcript_mrna_acc_type == 'CCDS':
                if best_transcript_mrna_acc_type in ['RefSeq', 'Ensembl']:
                    pass
                elif best_transcript_mrna_acc_type == 'CCDS':
                    if codon.aa_len > best_transcript_codon.aa_len:
                        best_transcript_codon = codon
        return best_transcript_codon

    def __get_codons(self, genomic):
        "Look up codons overlapping this genomic coordinate" 
        codon = None
        codons = []
        for row in CodonTblRow.query(self.conn, genomic):

            # If CodonTblRow.query returns any row,
            # then the position of "genomic" hits at least one
            # coding sequence region.

            codon = Codon(row)
            self.__set_codon_pos(genomic, codon)
            self.__set_ref_alt_codons(codon, genomic)
            self.__set_ref_alt_aas(codon, genomic)
            self.__set_sequence_ontology(codon, genomic)
       
            if codon.strand != genomic.strand:
                if codon.sequence_ontology in ['II', 'FI', 'CS']:
                    if codon.codon_pos == 0:
                        codon.codonNo -= 1
            codons.append(codon)

        # If multiple codon mappings are assigned to a genomic position,
        # the last one is chosen (by overwriting).

        #if len(codons) > 1:
        #    if self.__isAltSplice(codons):
        #        sys.stderr.write(variant_uid + "\t" + "alternatively_spliced_codon\n")
        #    else:
        #        sys.stderr.write(variant_uid + "\t" + "multiple_mappings\n")

        best_codon = self.__pick_representative_codon(codons)
        if self.pickOneTranscript == True:
            codons = [best_codon]

        return (codons, best_codon)

    def __isAltSplice(self, codons):
        "If at least 2 codons are split across an intron, case may be alt splice"
        altSplice = [x for x in codons.values() if x.size > 2]
        return len(altSplice) >= 1

    def __selectCodon(self, codons, genomic):
        "Arbitrarily select codon - prioritize missense"
        codonlist = []
        for codon in codons:
            alt_base_pos_in_codon = self.__get_codonsPos(genomic, codon)
            if self.__matches(genomic, codon):
                wtAA = self.__getRefAA(codon)
                pos = codon.codonNo
                altAA = self.__mutate(genomic, codon, alt_base_pos_in_codon)
                if wtAA != None and altAA != None:
                    if self.__isMissense(wtAA, altAA):
                        codonlist.append(codon)
        if len(codonlist) < 1:
            codonlist = codons
        # Commented by hcarter on 06/28/2009 - although less biased, the random selection can 
        # cause different mappings to be returned each time snvGet is run. Consistency is 
        # more reproducible and easier to deal with computationally.
        #codon = random.sample(codonlist,1)[0]
        #return codon
        return codonlist[0]
        
    def __getCodonPos(self, genomic, codon):
        """Determine position in the codon"""
        return [int(codon.pos1), int(codon.pos2), int(codon.pos3)].index(int(genomic.end))

    def __set_codon_pos(self, genomic, codon):
        """Determine position in the codon"""
        codon.codon_pos = [int(codon.pos1), int(codon.pos2), int(codon.pos3)].index(int(genomic.end))

    def __matches(self, genomic, codon):
        """Make sure that supplied reference base matches documented reference base"""
        return (genomic.strand == codon.strand and genomic.ref == codon.nucs[codon.codon_pos]) or (genomic.strand != codon.strand and rcDict.get(genomic.ref) == codon.nucs[codon.codon_pos])
    
    def __translate_codon_bases(self, codon_bases):
        bases = codon_bases.upper()
        if codonTable.has_key(bases):
            return codonTable.get(codon_bases.upper())
        else:
            return '?'

    def __getRefAA(self, codon):
        """Return amino acid before base substitution """
        return codonTable.get(codon.nucs)

    def __get_alt_codon(self, genomicVar, codon, alt_base_alt_base_pos_in_codon):
        alt_base = genomicVar.alt
        if genomicVar.strand != codon.strand:
            alt_base = rcDict.get(alt_base)
        alt_codon = codon.nucs[:alt_base_alt_base_pos_in_codon] + alt_base + codon.nucs[alt_base_alt_base_pos_in_codon+1:]
        return alt_codon

    def __mutate(self, genomicVar, codon, pos):
        """Return amino acid after base substitution """
        alt = genomicVar.alt
        if genomicVar.strand != codon.strand:
            alt = rcDict.get(alt)
        altCodon = codon.nucs[:pos] + alt + codon.nucs[pos+1:]
        return codonTable.get(altCodon)

    def __getMutationType(self):
        __getCdsMutationType(wtAA, altAA)

    def __getCdsMutationType(self, wtAA, altAA):
        """Returns coding sequence mutation type using Sequence Ontology terms """
        if wtAA == '*':
            if altAA == '*':
                mutationType = 'SY' #'synonymous_variant'
            else:
                mutationType = 'SL' #'stop_lost'
        elif altAA == '*':
            if wtAA == '*':
                mutationType = 'SY' #'synonymous_variant'
            else:
                mutationType = 'SG' #'stop_gained'
        else:
            if wtAA == altAA:
                mutationType = 'SY' #'synonymous_variant'
            else:
                mutationType = 'MS' #'missense_variant'
        return mutationType


class TranscriptMap(object):
    """
    Return a transcript and codon if it overlaps the genomic coordinate.
    """

    def __init__(self, dbconn):
        self.conn = dbconn

    def lookup(self, var):
        return CodonTblRow.transcriptquery(self.conn, var)

class DbRowBase(object):
    "base class used to convert a row of a query to an object"
    def __init__(self, row, description):
        for i in xrange(len(description)):
            if row != None:
                self.__dict__[description[i][0]] = row[i]
            else:
                self.__dict__[description[i][0]] = "NA"

    def get(self, colname):
        return self.__dict__.get(colname)

class CodonTblRow(DbRowBase):
    
    @staticmethod
    def query(conn, var):
        cursor=conn.cursor()
        cursor.execute('select t.aaLen as aa_len, t.RefseqT as refseq, t.EnsT as ens, t.CCDS as ccds, ct.* from (select * from CodonTable where chrom="' + var.chrom + '" and pos1=' + var.end + ' UNION select * from CodonTable where chrom="' + var.chrom + '" and pos2=' + var.end + ' UNION select * from CodonTable where chrom="' + var.chrom + '" and pos3=' + var.end + ') as ct, Transcript as t where ct.UID=t.UID')
        while True:
            row = cursor.fetchone()
            if row == None:
                break
            yield(CodonTblRow(row, cursor.description))

    @staticmethod
    def get_codon_bases(conn, uid, codon_no):
        cursor = conn.cursor()
        cursor.execute('select bases from CodonTable where UID=' + str(uid) + ' and Pos=' + str(codon_no))
        row = cursor.fetchone()
        if row == None:
            return None
        else:
            bases = row[0]
            return bases

    @staticmethod
    def transcriptquery(conn, var):
        cursor=conn.cursor()
        cursor.execute("Select UID from CodonTable where chrom=\""+var.chrom+"\" and pos1<="+var.end+" and pos3>="+var.end+";")
        row = cursor.fetchone()
        if row == None:
            transcript = None
        else:
            uid = row[0]
            cursor.execute('select RefseqT, CCDS, EnsT from Transcript where UID=' + str(uid))
            (dummy, ccds, refseq_t, dummy, ens_t, dummy, dummy) = cursor.fetchone()
            if refseq_t != None:
                transcript = refseq_t
            elif ccds != None:
                transcript = ccds
            elif ens_t != None:
                transcript = ens_t
        return transcript
        

class Codon(object):
    """Note: pos1, pos2, pos3 are 0-based"""
    def __init__(self, row):
        self.uid = row.get('UID')
        self.aa_len = row.get('aa_len')
        self.refseq = row.get('refseq')
        self.ens = row.get('ens')
        self.ccds = row.get('ccds')
        if self.refseq != None:
            self.mrnaAcc = self.refseq
            self.mrna_acc_type = 'RefSeq'
        elif self.ens != None:
            self.mrnaAcc = self.ens
            self.mrna_acc_type = 'Ensembl'
        elif self.ccds != None:
            self.mrnaAcc = self.ccds
            self.mrna_acc_type = 'CCDS'
        self.codonNo = row.get("Pos")
        self.chrom = row.get("chrom")
        self.pos1 = int(row.get("pos1"))
        self.pos2 = int(row.get("pos2"))
        self.pos3 = int(row.get("pos3"))
        self.nucs = row.get("bases")
        # Forward strand
        self.strand = "+"
        self.size = self.pos3-self.pos1
        # Reverse Strand
        if int(self.pos1) > int(self.pos3):
            self.strand = "-"
            self.size = self.pos1-self.pos3
        self.aa = codonTable[self.nucs]
   
    def __str__(self):
        return str(self.uid)+','+self.mrnaAcc+','+str(self.aa_len)+','+str(self.codonNo)+','+self.aa+','+self.chrom+','+self.strand+','+str(self.pos1)+','+str(self.pos2)+','+str(self.pos3)+','+self.nucs
