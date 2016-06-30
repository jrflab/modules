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

Created August 17, 2011

Contact: chasm-users@lists.johnshopkins.edu
'''

#External lib imports
import sys
import random
import string
import time

# Vars
AminoAcids = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

class FeatureFetcher(object):
    """
    Object to retrieve amino acid substitution features from the SNVBOX database
    """
    def __init__(self, dbconn, filename, rawValues=False, allowMissing=False):
        self.rawValues = rawValues
        self.allowMissing = allowMissing
        self.featureName = {}
        self.tableByFeature = {}
        self.featureMeans = {}
        self.featureRMSs = {}
        self.tableByName = {"AA_Features":AATable(),"Regional_Comp":RCTable(),"Local_Structure":LSTable(),"Genomic_MSA":GMSATable(),"Sam_MSA":SMSATable(),"Exon_Features":ExonTable(),"Uniprot_features":UniprotTable(),"Triplet_Features":TripletTable(),"Triplet_Sub_Features":TripletSubTable()} # enable retrieval of table object by table name
        self.orderedList = ["AA_Features","Genomic_MSA","Sam_MSA","Exon_Features","Local_Structure","Regional_Comp","Triplet_Features","Triplet_Sub_Features","Uniprot_features"]
        self.orderedFeatureList = []
        self.orderedFeatureNames = []
        self.__getSnvBoxFeatures(dbconn)
        self.__populateFromFeatureList(filename)
        self.__buildOrderedList()

    def __getSnvBoxFeatures(self, dbconn):
        """Store SNVbox feature mappings (SNVBox Table and column of table containing feature values)"""
        for row in FeatureTblRow.query(dbconn):
            self.featureName[row.get("Feature")] = row.get("Col_Name")
            self.tableByFeature[row.get("Feature")] = row.get("Table_Name")
            self.featureMeans[row.get("Feature")] = float(row.get("Mean"))
            self.featureRMSs[row.get("Feature")] = float(row.get("RMS"))
            
    def __populateFromFeatureList(self, filename):
        """
        Populate feature table objects with the features the user specified
        - Only lookup requested features
        """
        for feature in FeatureFileReader(filename):
            tableName = self.tableByFeature.get(feature)
            if tableName != None: 
                table = self.tableByName.get(tableName)
                table.featurelist.append(self.featureName.get(feature))
                table.featureMean[feature] = self.featureMeans.get(feature)
                table.featureRMS[feature] = self.featureRMSs.get(feature)
                table.featurenames.append(feature)
            else:
                raise Exception(feature + ": No such feature in SNVBox\n")

    def __buildOrderedList(self):
        """Assign feature order - must be consistent for all arffs generated"""
        for tableName in self.orderedList:
            self.orderedFeatureNames += self.tableByName.get(tableName).featurenames
            self.orderedFeatureList += self.tableByName.get(tableName).featurelist

    def getFeatures(self, dbconn, var):
        """
        Retrieve features for variant object 
        -if rawValues=True return unscaled feature values 
        -if allowMissing=True, return None for missing values
        """
        if var.alt_aa not in AminoAcids or var.ref_aa not in AminoAcids:
            return None
        q = TranscriptFeatures.query(dbconn, var.acc)
        var.UID = q.get("UID")
        if var.UID != None:
            var.UID = str(var.UID)
            tsID = q.get("tsID")
            if self.__isTsVerMismatch(var.acc, tsID):
                sys.stderr.write("Warning: Refseq transcript version number (" + var.acc + ") does not match Refseq version in database (" + tsID + ")\n")
            feats = []
            featNames = []
            for tableName in self.orderedList:
                if len(self.tableByName.get(tableName).featurelist) > 0:
                    feats += self.tableByName.get(tableName).getFeatures(dbconn, var)
                    featNames += self.tableByName.get(tableName).featurenames
            #print '--- before scaling and missing value filling ---'
            #print var.UID, feats
            if not self.rawValues:
                feats = self.__scaleValues(feats, featNames)
                #print '!!! after scaling !!!'
                #print var.UID, feats
            if not self.allowMissing and None in feats:
                feats = self.__fillMissingValues(feats)
                #print '+++ after filling missing values +++'
                #print var.UID, feats
        else:
            if var.genomicstatus == 1:
                sys.stderr.write(var.uid+" "+var.chrom+":"+var.start+"-"+var.end+" "+var.strand+" " +var.ref_aa+">"+var.alt_aa+" mapped to "+var.acc+" "+var.aaSubst+" was not found in SNVbox database\n")
            else:
                sys.stderr.write(var.uid + " " + var.acc + " " + var.ref_aa+var.aa_pos+var.alt_aa + " not found in SNVbox database\n")
            feats = None
        return feats
    
    def __isTsVerMismatch(self, varID, dbID):
        """Check for refseq version mismatch """
        mismatch = False
        if varID.find(".") != -1 and not varID.startswith("CCDS"):
            ver = varID.split(".")[-1]
            dbver = dbID.split(".")[-1]
            if ver != dbver:
                mismatch = True
        return mismatch
            
    def __scaleValues(self, feats, featNames):
        """Scale features by subtracting mean and dividing by RMS"""
        for i in xrange(len(feats)):
            featMean = self.featureMeans.get(featNames[i])
            featRMS = self.featureRMSs.get(featNames[i])
            if feats[i] != None:
                if featRMS == 0:
                    feats[i] = 0.0 
                else:
                    feats[i] = str((float(feats[i])-featMean)/featRMS)
        return feats

    def __fillMissingValues(self, feats):
        """Fill missing values - can fill with mean raw value or mean after scaling ( =0 )"""
        for i in xrange(len(feats)):
            if feats[i] == None:
                if self.rawValues:
                    feats[i] = str(self.featureMeans.get(self.orderedFeatureList[i]))
                else:
                    feats[i] = "0.0" # Assumes mean value post-scaling is 0
        return feats


class FeatureFileReader(object):
    '''
    Object to read featurelist file
    '''
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
            return line[:-1]


class DbTableBase(object):
    """Base class for storing feature information retrieved from SNVBOX tables"""
    def __init__(self):
        self.featurelist = []
        self.featureValues = []
        self.featurenames = []
        self.featureMean = {}
        self.featureRMS = {}

class AATable(DbTableBase):
    
    def getFeatures(self, conn, var):
        self.featureValues = []
        row = AAFeatures.query(conn, self.featurelist, var.ref_aa, var.alt_aa)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues
            
class RCTable(DbTableBase):

    def getFeatures(self, conn, var):
        self.featureValues = []
        row = RCFeatures.query(conn, self.featurelist, var.UID, var.aa_pos)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class LSTable(DbTableBase):

    def getFeatures(self, conn, var):
        self.featureValues = []
        row = LSFeatures.query(conn, self.featurelist, var.UID, var.aa_pos)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class GMSATable(DbTableBase):
    
    def getFeatures(self, conn, var):
        self.featureValues = []
        row = GMSAFeatures.query(conn, self.featurelist, var.UID, var.aa_pos, var.alt_aa)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class SMSATable(DbTableBase):

    def getFeatures(self, conn, var):
        self.featureValues = []
        row = SMSAFeatures.query(conn, self.featurelist, var.UID, var.aa_pos, var.alt_aa)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class ExonTable(DbTableBase):

    def getFeatures(self, conn, var):
        self.featureValues = []
        row = ExonFeatures.query(conn, self.featurelist, var.UID, var.aa_pos)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class UniprotTable(DbTableBase):

    def getFeatures(self, conn, var):
        self.featureValues = []
        row = UniprotFeatures.query(conn, self.featurelist, var.UID, var.aa_pos)
        for feature in self.featurelist:
            self.featureValues.append(row.get(feature))
        return self.featureValues

class Triplets:
    uid = None
    pos = None
    triplets = None

    @staticmethod
    def findTriplets(conn, var):
        cursor = conn.cursor()

        Triplets.uid = var.UID
        Triplets.pos = var.aa_pos

        # Gets the three triplets (having the mutation amino acid in the first, second, and third position).
        first_triplet = None
        second_triplet = None
        third_triplet = None
        pos = int(var.aa_pos)
        pos_m2 = pos - 2
        pos_m1 = pos - 1
        pos_p1 = pos + 1
        pos_p2 = pos + 2

        cursor.execute('select Pos, WildType from Regional_Comp where UID=' + var.UID + ' and Pos>=' + str(pos_m2) + ' and Pos<=' + str(pos_p2))
        results = cursor.fetchall()
        pos_to_aa = {}
        for result in results:
            (pos, aa) = result
            pos_to_aa[pos] = aa
        
        if pos_to_aa.has_key(pos_m2) == False or pos_to_aa.has_key(pos_m1) == False or pos_to_aa.has_key(pos) == False:
            triplet_3_wt = None
            triplet_3_mut = None
        else:
            triplet_3_wt = pos_to_aa[pos_m2] + pos_to_aa[pos_m1] + var.ref_aa
            triplet_3_mut = pos_to_aa[pos_m2] + pos_to_aa[pos_m1] + var.alt_aa

        if pos_to_aa.has_key(pos_m1) == False or pos_to_aa.has_key(pos) == False or pos_to_aa.has_key(pos_p1) == False:
            triplet_2_wt = None
            triplet_2_mut = None
        else:
            triplet_2_wt = pos_to_aa[pos_m1] + var.ref_aa + pos_to_aa[pos_p1]
            triplet_2_mut = pos_to_aa[pos_m1] + var.alt_aa + pos_to_aa[pos_p1]

        if pos_to_aa.has_key(pos) == False or pos_to_aa.has_key(pos_p1) == False or pos_to_aa.has_key(pos_p2) == False:
            triplet_1_wt = None
            triplet_1_mut = None
        else:
            triplet_1_wt = var.ref_aa + pos_to_aa[pos_p1] + pos_to_aa[pos_p2]
            triplet_1_mut = var.alt_aa + pos_to_aa[pos_p1] + pos_to_aa[pos_p2]

        Triplets.triplets = [triplet_1_wt, triplet_1_mut, triplet_2_wt, triplet_2_mut, triplet_3_wt, triplet_3_mut]

class TripletTable(DbTableBase):

    def getFeatures(self, conn, var):
        if Triplets.uid != var.UID or Triplets.pos != var.aa_pos:
            Triplets.findTriplets(conn, var)
        [triplet_1_wt, triplet_1_mut, triplet_2_wt, triplet_2_mut, triplet_3_wt, triplet_3_mut] = Triplets.triplets

        #cursor = conn.cursor()
        #
        # Gets the three triplets (having the mutation amino acid in the first, second, and third position).
        #first_triplet = None
        #second_triplet = None
        #third_triplet = None
        #pos = int(var.aa_pos)
        #pos_m2 = pos - 2
        #pos_m1 = pos - 1
        #pos_p1 = pos + 1
        #pos_p2 = pos + 2
        #
        #cursor.execute('select Pos, WildType from Regional_Comp where UID=' + var.UID + ' and Pos>=' + str(pos_m2) + ' and Pos<=' + str(pos_p2))
        #results = cursor.fetchall()
        #pos_to_aa = {}
        #for result in results:
        #    (pos, aa) = result
        #    pos_to_aa[pos] = aa
        #
        #if pos_to_aa.has_key(pos_m2) == False or pos_to_aa.has_key(pos) == False:
        #    triplet_3_wt = None
        #    triplet_3_mut = None
        #else:
        #    triplet_3_wt = pos_to_aa[pos_m2] + pos_to_aa[pos_m1] + var.ref_aa
        #    triplet_3_mut = pos_to_aa[pos_m2] + pos_to_aa[pos_m1] + var.alt_aa
        #
        #if pos_to_aa.has_key(pos_m1) == False or pos_to_aa.has_key(pos_p1) == False:
        #    triplet_2_wt = None
        #    triplet_2_mut = None
        #else:
        #    triplet_2_wt = pos_to_aa[pos_m1] + var.ref_aa + pos_to_aa[pos_p1]
        #    triplet_2_mut = pos_to_aa[pos_m1] + var.alt_aa + pos_to_aa[pos_p1]
        #
        #if pos_to_aa.has_key(pos) == False or pos_to_aa.has_key(pos_p2) == False:
        #    triplet_1_wt = None
        #    triplet_1_mut = None
        #else:
        #    triplet_1_wt = var.ref_aa + pos_to_aa[pos_p1] + pos_to_aa[pos_p2]
        #    triplet_1_mut = var.alt_aa + pos_to_aa[pos_p1] + pos_to_aa[pos_p2]

        # The order of self.featurelist is ['FirstProbMut', 'FirstProbWild', 'SecondProbMut', 'SecondProbWild', 'ThirdProbMut', 'ThirdProbWild']
        self.featureValues = []
        if triplet_1_mut == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_1_mut)
            if prob != None:
                prob = prob[0]
            self.featureValues.append(prob)
        if triplet_1_wt == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_1_wt)
            if prob != None:
                prob = prob[0]
            self.featureValues.append(prob)
        if triplet_2_mut == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_2_mut)
            if prob != None:
                prob = prob[1]
            self.featureValues.append(prob)
        if triplet_2_wt == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_2_wt)
            if prob != None:
                prob = prob[1]
            self.featureValues.append(prob)
        if triplet_3_mut == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_3_mut)
            if prob != None:
                prob = prob[2]
            self.featureValues.append(prob)
        if triplet_3_wt == None:
            self.featureValues.append(None)
        else:
            prob = TripletFeatures.query(conn, triplet_3_wt)
            if prob != None:
                prob = prob[2]
            self.featureValues.append(prob)

        #for i in xrange(len(self.featurelist)):
        #    print 'i=',i,', feature=',self.featurelist[i]
        #    if self.featurenames[i] in ["AATripletFirstProbWild","AATripletSecondProbWild","AATripletThirdProbWild"]:
        #        row1 = TripletFeatures.query(conn, self.featurelist[i], var.UID, var.aa_pos, var.ref_aa)
        #        self.featureValues.append(row1.get(self.featurelist[i]))
        #    elif self.featurenames[i] in ["AATripletFirstProbMut","AATripletSecondProbMut","AATripletThirdProbMut"]:
        #        row2 = TripletFeatures.query(conn, self.featurelist[i], var.UID, var.aa_pos, var.alt_aa)
        #        self.featureValues.append(row2.get(self.featurelist[i]))
        return self.featureValues

class TripletSubTable(DbTableBase): 

    def getFeatures(self, conn, var):
        if Triplets.uid != var.UID or Triplets.pos != var.aa_pos:
            Triplets.findTriplets(conn, var)
        [triplet_1_wt, triplet_1_mut, triplet_2_wt, triplet_2_mut, triplet_3_wt, triplet_3_mut] = Triplets.triplets

        self.featureValues = []
        if triplet_1_wt == None:
            self.featureValues.append(None)
        else:
            [firstdiffprob, seconddiffprob, thirddiffprob] = TripletSubFeatures.query(conn, triplet_1_wt[1:3], var.ref_aa, var.alt_aa)
            self.featureValues.append(firstdiffprob)
        if triplet_2_wt == None:
            self.featureValues.append(None)
        else:
            [firstdiffprob, seconddiffprob, thirddiffprob] = TripletSubFeatures.query(conn, triplet_2_wt[0] + triplet_2_wt[2], var.ref_aa, var.alt_aa)
            self.featureValues.append(seconddiffprob)
        if triplet_3_wt == None:
            self.featureValues.append(None)
        else:
            [firstdiffprob, seconddiffprob, thirddiffprob] = TripletSubFeatures.query(conn, triplet_3_wt[:2], var.ref_aa, var.alt_aa)
            self.featureValues.append(thirddiffprob)
        #for feature in self.featurelist:
        #    row = TripletSubFeatures.query(conn, feature, var.UID, var.aa_pos, var.ref_aa, var.alt_aa)
        #    self.featureValues.append(row.get(feature))
        return self.featureValues

class DbRowBase(object):
    "base class used to convert a row of a query to an object"
    def __init__(self):
        pass

    def __init__(self, row, description):
        for i in xrange(len(description)):
            if row != None:
                self.__dict__[description[i][0]] = row[i] if row[i] != None else 0.0
            else:
                self.__dict__[description[i][0]] = None

    def get(self, colname):
        return self.__dict__.get(colname)
            
    
class FeatureTblRow(DbRowBase):

    @staticmethod
    def query(conn):
        cursor=conn.cursor()
        cursor.execute("Select * from Features;")
        while True:
            row = cursor.fetchone()
            if row == None:
                break
            yield(FeatureTblRow(row, cursor.description))
        

class TranscriptFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, acc):
        query = TranscriptFeatures.buildQuery(acc)
        cursor = conn.cursor()
        cursor.execute(query)
        row = cursor.fetchone()
        feats = TranscriptFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                pass
                # Commented by HC (11/18/2012)
                # Can have multiple transcripts with the same CDS so multiple
                # query matches should be allowed.
                # raise Exception("Unexpected duplicate hit in database\n")
            cursor.close()
        return feats

    @staticmethod
    def buildQuery(acc):
        if acc.find(".") != -1:
            core = acc.split(".")[0]
            ver = acc.split(".")[1]
        else:
            core = acc
            ver = None
        if acc.startswith("NP"):
            return "Select UID, RefseqP as tsID from Transcript where RefseqP like \"" + core + ".%\";"
        elif acc.startswith("NM"):
            return "Select UID, RefseqT as tsID from Transcript where RefseqT like \"" + core + ".%\";"
        elif acc.startswith("ENSP"):
            return "Select UID, EnsP as tsID from Transcript where EnsP=\""+ core +"\";"
        elif acc.startswith("ENST"):
            return "Select UID, EnsT as tsID from Transcript where EnsT=\""+ core +"\";"
        elif acc.startswith("CCDS"):
            return "Select UID, CCDS as tsID from Transcript where CCDS like \""+ core +".%\";"
        else:
            raise Exception(acc + " not recognized. Accession may not supported by SNVbox\n")

class AAFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, wt, alt):
        cursor = conn.cursor()
        cursor.execute("Select " + ",".join(features) + " from AA_Features where WildType=\"" + wt+ "\" and Mut = \"" + alt + "\";")
        row = cursor.fetchone()
        feats = AAFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                sys.stderr.write("AAFeatures:Duplicate hit:"+stmt+"\n")
                #raise Exception("Unexpected duplicate hit in database\n")
            cursor.close()
        return feats

class RCFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos):
        cursor = conn.cursor()
        cursor.execute("Select " + ",".join(features) + " from Regional_Comp where UID = \"" + UID + "\" and Pos= " + pos+ ";")
        row = cursor.fetchone()
        feats = RCFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                sys.stderr.write("RCFeatures:Duplicate hit:"+stmt+"\n")
                #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
            cursor.close()
        return feats

class LSFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos):
        cursor = conn.cursor()
        cursor.execute('select LocalStructureId from CodonToLocalStructure where UID=' + UID + ' and Pos=' + pos)
        result = cursor.fetchone()
        if result != None:
            localstructureid = result[0]
            cursor.execute('select ' + ','.join(features) + ' from Local_Structure where id=' + str(localstructureid))
            row = cursor.fetchone()
            feats = LSFeatures(row,cursor.description)
            if row != None:
                row = cursor.fetchone()
                if row != None:
                    sys.stderr.write("LSFeatures:Duplicate hit:"+stmt+"\n")
                    #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
                cursor.close()
        else:
            feats = LSFeatures(None, cursor.description)
        return feats

class GMSAFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos,  alt):
        q_features = [x for x in features if x != "PHC"]
        if "PHC" in features:
            q_features.append("PHC_" + alt + " as PHC")                    
        cursor = conn.cursor()
        cursor.execute("Select " + ",".join(q_features) + " from Genomic_MSA where UID = \"" + UID + "\" and Pos= " + pos+ ";")
        row = cursor.fetchone()
        feats = GMSAFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                sys.stderr.write("GMSAFeatures:Duplicate hit:"+stmt+"\n")
                #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
            cursor.close()
        return feats

class SMSAFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos, alt):
        q_features = [x for x in features if x != "PHC"]
        if "PHC" in features:
            q_features.append("PHC_" + alt + " as PHC")
        cursor = conn.cursor()
        cursor.execute('select SamMSAId from CodonToSamMSA where UID=' + UID + ' and Pos=' + pos)
        result = cursor.fetchone()
        if result != None:
            sammsaid = result[0]
            cursor.execute('select ' + ','.join(q_features) + ' from Sam_MSA where id=' + str(sammsaid))
            row = cursor.fetchone()
            feats = SMSAFeatures(row,cursor.description)
            if row != None:
                row = cursor.fetchone()
                if row != None:
                    sys.stderr.write("SMSAFeatures:Duplicate hit:"+stmt+"\n")
                    #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
                cursor.close()
        else:
            feats = SMSAFeatures(None, cursor.description)
        return feats
        
class ExonFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos):
        q_features = ["e." + x for x in features]
        cursor = conn.cursor()
        stmt = "Select " + ",".join(q_features) + " from Transcript_Exon t, Exon_Features e where t.tStart <= "+str((int(pos)- 1)*3)+" and t.tend >= "+str((int(pos) - 1)*3)+" and t.UID=e.UID and t.Exon=e.Exon and t.UID="+UID+";"
        cursor.execute(stmt)
        row = cursor.fetchone()
        feats = ExonFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                sys.stderr.write("ExonFeatures:Duplicate hit:"+stmt+"\n")
            cursor.close()
        return feats

class UniprotFeatures(DbRowBase):
    ""

    @staticmethod
    def query(conn, features, UID, pos):
        q_features = ["u." + x for x in features]
        cursor = conn.cursor()
        cursor.execute("Select " + ",".join(q_features) + " from Uniprot_Xref x left join  Uniprot_features u on u.Acc=x.Uniprot and u.Pos=x.UniprotPos where x.UID="+ UID +" and x.Pos="+ pos +";")
        row = cursor.fetchone()
        feats = UniprotFeatures(row,cursor.description)
        if row != None:
            row = cursor.fetchone()
            if row != None:
                sys.stderr.write("UniprotFeatures:Duplicate hit:"+stmt+"\n")
                #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
            cursor.close()
        return feats
        
class TripletFeatures(DbRowBase):
    ""

    tripletFeaturesData = None

    @staticmethod
    def query(conn, triplet): 
        cursor = conn.cursor()

        # Loads Triplet_Features table into memory.
        if TripletFeatures.tripletFeaturesData == None:
            TripletFeatures.tripletFeaturesData = {}
            cursor.execute('select * from Triplet_Features')
            results = cursor.fetchall()
            for result in results:
                (ttt, firstprob, secondprob, thirdprob) = result
                TripletFeatures.tripletFeaturesData[ttt] = [firstprob, secondprob, thirdprob]

        if TripletFeatures.tripletFeaturesData.has_key(triplet):
            return TripletFeatures.tripletFeaturesData[triplet]
        else:
            return None

        # Col_names in Feature table are FirstProbMut, FirstProbWild, but actual column name in Triplet_Features table is FirstProb -> need to drop Mut or Wild
        #if feature.endswith("Mut"):
        #    feature = feature[:-3]
        #    suffix = "Mut"
        #else:
        #    feature = feature[:-4]
        #    suffix = "Wild"
        #mkTriplet = TripletFeatures.formatQuery(feature, UID, pos, res)
        #cursor.execute("Select t." + feature + " as " + feature + suffix + " from Triplet_Features t, Regional_Comp r1, Regional_Comp r2 where t.Triplet = " + mkTriplet+ ";")
        #print "Select t." + feature + " as " + feature + suffix + " from Triplet_Features t, Regional_Comp r1, Regional_Comp r2 where t.Triplet = " + mkTriplet+ ";"
        #row = cursor.fetchone()
        #feats = TripletFeatures(row,cursor.description)
        #if row != None:            
        #    row = cursor.fetchone()
        #    if row != None:
        #        sys.stderr.write("TripletFeatures:Duplicate hit:"+stmt+"\n")
        #        #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
        #    cursor.close()
        #return feats

    @staticmethod
    def formatQuery(feature, UID, pos, res):
        if feature == "FirstProb":
            mkTriplet = "CONCAT(\""+res+"\", r1.WildType, r2.WildType) and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"+1 and r2.Pos="+pos+"+2"
        elif feature == "SecondProb":
            mkTriplet = "CONCAT(r1.WildType, \""+res+"\", r2.WildType) and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"+-1 and r2.Pos="+pos+"+1"
        elif feature == "ThirdProb":
            mkTriplet = "CONCAT(r1.WildType, r2.WildType, \""+res+"\") and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"-2 and r2.Pos="+pos+"-1"
        return mkTriplet

class TripletSubFeatures(DbRowBase):
    ""
    
    tripletSubFeatures = None

    @staticmethod
    def query(conn, two, wt, mut):
        cursor = conn.cursor()
        
        # Loads the whole TripletSubFeatures table.
        if TripletSubFeatures.tripletSubFeatures == None:
            TripletSubFeatures.tripletSubFeatures = {}
            cursor.execute('select * from Triplet_Sub_Features')
            results = cursor.fetchall()
            for result in results:
                (twobases, wildtype, mutant, firstdiffprob, seconddiffprob, thirddiffprob) = result
                if TripletSubFeatures.tripletSubFeatures.has_key(twobases) == False:
                    TripletSubFeatures.tripletSubFeatures[twobases] = {}
                if TripletSubFeatures.tripletSubFeatures[twobases].has_key(wildtype) == False:
                    TripletSubFeatures.tripletSubFeatures[twobases][wildtype] = {}
                TripletSubFeatures.tripletSubFeatures[twobases][wildtype][mutant] = [firstdiffprob, seconddiffprob, thirddiffprob]

        if TripletSubFeatures.tripletSubFeatures.has_key(two):
            if TripletSubFeatures.tripletSubFeatures[two].has_key(wt):
                if TripletSubFeatures.tripletSubFeatures[two][wt].has_key(mut):
                    return TripletSubFeatures.tripletSubFeatures[two][wt][mut]

        return [None, None, None]

        #mkTwoBases = TripletSubFeatures.formatQuery(feature, UID, pos)
        #cursor.execute("Select t." + feature + " from Triplet_Sub_Features t, Regional_Comp r1, Regional_Comp r2 where t.TwoBases="+ mkTwoBases+" and t.WildType=\""+wt+"\" and t.Mutant=\""+alt+"\";")
        #row = cursor.fetchone()
        #feats = TripletSubFeatures(row,cursor.description)
        #if row != None:            
        #    row = cursor.fetchone()
        #    if row != None:
        #        sys.stderr.write("TripletSubFeatures:Duplicate hit:"+stmt+"\n")
        #        #raise Exception("Unexpected duplicate hit in database for "+UID+" "+pos+"\n")
        #    cursor.close()
        #return feats

    @staticmethod
    def formatQuery(feature, UID, pos):
        if feature == "FirstDiffProb":
            mkTwoBases = "CONCAT(r1.WildType, r2.WildType) and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"+1 and r2.Pos="+pos+"+2"
        elif feature == "SecondDiffProb":
            mkTwoBases = "CONCAT(r1.WildType, r2.WildType) and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"+-1 and r2.Pos="+pos+"+1"
        elif feature == "ThirdDiffProb":
            mkTwoBases = "CONCAT(r1.WildType, r2.WildType) and r1.UID="+UID+" and r2.UID="+UID+" and r1.Pos="+pos+"-2 and r2.Pos="+pos+"-1"
        return mkTwoBases


if __name__ == "__main__":
    pass




