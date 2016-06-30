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

'''
Config.py

A class for reading simple key=value configuration file strings

'''

#Python imports
import os

class Config:
    '''
    This class implements a wrapper around configuration files. It attempts to
    read two files from the directory specified by the HOME environment variable.
    If home is not defined, then the instantiation of this object fails.
    During initialization, it reads .hg.conf and .lssnp.conf and loads all key/
    value pairs from these files into a dictionary called dict_. Keys in the
    .lssnp.conf will take precedence over identical keys in .hg.conf.
    '''
    
    #The configuration attribute
    dict_ = {}
    
    def __init__(self):
        '''
        Initialize the instance by reading in the configuration from the user's
        home directory. This class uses the values in .hg.conf first, followed
        by any values that exist in .lssnp.conf. The latter will override the 
        former if both files define the same key.
        '''
        #Get a path to the configuration files
        configPath = self._getConfigPath()
        self._loadConfigs(os.path.join(configPath,"snv_box.conf"))
    
    def get(self, key):
        '''
        Retrieves a value from the configuration dictionary
        '''
        return self.dict_.get(key)


    def getWDefault(self, key, default):
        '''
        Retrieves a value from the configuration dictionary, returning a default
        value if the key is not found.      
        '''
        return self.dict_.get(key, default)
    
    
    def _loadConfigs(self, path):
        '''
        This function loads into the config_ dictionary all of the key/value
        pairs found in the specified file.        
        '''
        try:
            #Open the configuration file
            f = open(path, 'r')
            
            try:
                #For each line in the file, see if we have an '=' and parse a
                # key/value pair into the config_ dictionary
                for line in f:
                    if line.find('=') >= 0:
                        key, value = line.split('=', 2)
                        self.dict_[key] = value.strip()
            finally:
                #Close the file
                f.close()
                    
        except IOError:
            msg = "Config:_loadConfigs - IOError: Unable to open configuration file "+path
            raise Exception(msg)


    def _getConfigPath(self):
        '''
        Retrieves the configuration path from the environment. This function 
        will PREFER the value of the HOME variable, if defined. If this variable
        is not defined, it will look for the USERNAME and possibly USERDOMAIN
        variables and create a /home/USERNAME.USERDOMAIN path. If the CYGWIN
        variable is defined, then we will assume that any home directory is
        within the CYGWIN root directory. We attempt to locate this in C:\cygwin.       
        '''
        
        # Reads from the source directory of the code
        home = os.path.join(os.environ['SNVBOXDIR'])
        return home + os.sep
    
    
if __name__ == "__main__":
    # Quick test
    config = Config()
    #print config.get("dataDir")

