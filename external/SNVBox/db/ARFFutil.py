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

Created August 16, 2011

Contact: chasm-users@lists.johnshopkins.edu
'''


# Vars
arff_dict = {}

# Classes

class ARFFwriter(object):
    """Object to write arff files in format used by parf """

    def __init__(self, arffname):
        self.arffhandle = open(arffname,"w")
        self.arffhandle.write("@relation headerfile\n\n")
        self.arffhandle.write("@attribute UID string\n")
        self.arffhandle.write("@attribute ID string\n")

    def addFeature(self, feature):
        self.arffhandle.write("@attribute " + feature + " numeric\n")

    def addClass(self, classlabels, weights=None):
        # Must be called explicity for arff to include a class definition in the header
        weights = ["(1)" for x in classlabels] if weights == None else ["("+str(x)+")" for x in weights]
        if len(weights) != len(classlabels):
            raise Exception("Number of weights must match number of class labels")
        classlabels = zip(classlabels, weights)
        labels = [" ".join(x) for x in classlabels]
        self.arffhandle.write("@attribute CLASS { " + ",".join(labels) + " }\n")

    def dataStart(self):
        self.arffhandle.write("\n@data\n\n")

    def addData(self, uid, name, values, classlabel=None):
        if values != None:
            if classlabel != None:
                outstr = self.__manageLength(str(uid) + " " + name + " " + " ".join(map(str,values)) + " " + classlabel + "\n")
                self.arffhandle.write(outstr)
            else:
                outstr = self.__manageLength(str(uid) + " " + name + " " + " ".join(map(str,values)) + "\n")
                self.arffhandle.write(outstr)

    def __manageLength(self, text):
        """PARF cannot take text strings longer than 1024 characters. Long strings must be split over multiple lines using & """
        if len(text) > 1024:
            index = 1000
            while index < len(text):
                while text[index] != " ":
                    index -= 1
                text = text[0:index] + " &\n" + text[index:]
                index += 1000
        return text

    def closeArff(self):
        self.arffhandle.close()



if __name__ == "__main__":
    
    # Test code
    testArff = ARFFwriter("arff.test")
    testArff.addFeature("testFeat")
    testArff.addClass(["testClass1","testClass2"])
    testArff.dataStart()
    testArff.addData("testID","testName",["testValue"],"testClass1")
    testArff.closeArff()
