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

Created August 18, 2011

Contact: chasm-users@lists.johnshopkins.edu
'''

#External lib imports
import MySQLdb

#Local imports
from db import Config

class DBUtil:
    '''
    Object to manage MySQLdb connections 
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.config_ = Config.Config()
        self.chasmDB_ = self.config_.getWDefault("chasmDB","SNVBox")
        self.conn = self._connect(self.chasmDB_)

    def close(self):
        self._release(self.conn)

    def _release(self, conn):
        '''
        Release a connection to the database. 
        '''
        conn.close();

    def _connect(self, dbname):
        '''
        Establish a connection to the database. This functions uses credentials
        defined in the configuration and the database name provided. It returns
        a valid connection or throws an exception.
        '''
        #Get configuration instance
        conn = None

        try:
            #Load all db connect data from configuration
            dbuser = self.config_.get("db.user")
            dbpassword = self.config_.get("db.password")
            dbhost = self.config_.get("db.host")

            dbsocket = self.config_.getWDefault("db.unix_socket", "")
            dbport = int(self.config_.getWDefault("db.port","3306"))

            if dbhost == "localhost" and dbsocket!="":
                # Connect using unix socket
                conn = MySQLdb.connect(host=dbhost,
                                   user=dbuser,
                                   passwd=dbpassword,
                                   db=dbname,
                                   unix_socket=dbsocket)
            else:
                # Connect using TCPIP
                conn = MySQLdb.connect(host=dbhost,
                                   user=dbuser,
                                   passwd=dbpassword,
                                   db=dbname,
                                   port=dbport)

        except MySQLdb.Error,e:
            print """Error %d: %s""" % (e.args[0], e.args[1])
            print "Unable to connect to database, please check snv_box.conf settings."
            print "Currently SnvGet is using the following settings:"
            print "\tdb.user     = " + dbuser
            print "\tdb.password = " + dbpassword
            print "\tdb.host     = " + dbhost
            print "\tdb.port     = " + str(dbport)
            print "\tdb.unix_socket = " + dbsocket
            print ""
            print "By default, SnvGet connects using TCP/IP unless host is set to"
            print "\"localhost\" and db.unix_socket is defined. If neither db.port nor"
            print "db.unix_socket is defined, it will try connect using port 3306."
        except Exception, inst:
            msg = "DBUtil:_connect - Unknown error caught connecting to database"
            raise Exception(str(type(inst)), str(inst), msg)

        #Return the connection
        return conn

