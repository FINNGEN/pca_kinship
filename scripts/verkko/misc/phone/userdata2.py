"""
Module for analysing mobile phone network metadata.

NOTE: This file is used for mobile phone data set 2. For set 1, use
userdata.py instead.

This file works with the new data (set2) received in November 2009.

Available fields:

  'age'            Customer age. 0 = missing/unknown
  'gender'         Customer gender. The value is either UserData.male,
                   UserData.female or UserData.unknown (= 0). If a
                   non-postpaid user is listed as male and has either
                   no valid ZIP or no valid age, the gender will be
                   marked as unknown. This is done because male
                   appears to be the default gender.
  'ZIP'            Customer ZIP code. Note that this is save as integer
                   for performance reasons, use formatting '%05d' to 
                   zero-pad to proper length.
  'validZIP'       True if ZIP is valid. ZIP is valid if it is included in
                   the list of valid ZIP codes in the country in question.
                   Corrupted ZIP codes (non-numeral characters, less than 
                   five digits, etc.) are always invalid.
  'usertype'       Either prepaid or postpaid. The value is either
                   UserData.prepaid, UserData.postpaid or 
                   UserData.unknown (= 0).
  'contractID'     Customer's contract ID. Several customers can have the
                   same contract ID.
  'activationDate' The date when the phone number was first activated.
                   If empty, the activation date is unknown.
  'disconnectDate' The date when the phone number was first disconnected.
                   If empty, the customer is still active.
  'latitude'       Latitude of the most most common position, 0 if missing.
  'longitude'      Longitude of the most most common position, 0 if missing.


Usage:

>>> # Read user data either from a text file or from a previously
>>> # saved npy-file. Reading from npy-file is _much_ faster!
>>> from phone.userdata2 import UserData
>>> udata = UserData('userdata.txt')
>>> udata.save('userdata.npy')
>>> del udata
>>> udata = UserData('userdata.npy')

>>> # The original ASCII file name is also stored:
>>> udata.originalFileName()
'userdata.txt'

>>> # Take a look at the data of a single user.
>>> udata[1127]
(52, 1, 40297L, True, 2)
>>> # The previous line means that user 1127 is 52 year old female
>>> # (the second '1'), with valid ZIP 40297 (validZIP = True), and
>>> # is a postpaid user (the last '2').

>>> # User data should be compared to internal constants. For example,
>>> # 'gender' is one of UserData.female, UserData.male or
>>> # UserData.unknown:
>>> udata[1127]['gender'] == UserData.female
True
>>> udata[1127]['gender'] == UserData.male
False
>>> udata[1127]['gender'] == UserData.unknown
False

>>> # 'usertype' is one of UserData.postpaid, UserData.prepaid or
>>> # UserData.unknown:
>>> udata[1127]['usertype'] == UserData.postpaid
True
>>> udata[1127]['usertype'] == UserData.prepaid
False
>>> udata[1127]['usertype'] == UserData.unknown
False

>>> # Test if ZIP code is valid (see `help UserData` for condition of
>>> # validity). Note also that ZIP is saved as integer; use
>>> # formatting '%05d' to print it properly.
>>> if udata[1127]['validZIP']:
...     print "ZIP %05d is valid!" % udata[1127]['ZIP']
...
ZIP 40297 is valid!

>>> # Type UserData.unknown evaluates to False. Therefore it is
>>> # possible to go through valid users like this:
>>> N_valid_users = 0
>>> for user in udata:
...     if (user['gender'] and user['usertype'] and user['validZIP']):
...         N_valid_users += 1
...
>>> print "Only %d valid users out of %d!" % (N_valid_users, len(udata))
Only 4162538 valid users out of 5343749!


"""

import os
import numpy as np

def getValidZIPs():
    """Return a set of all valid ZIP codes.
    
    Return
    ------
    zips : set of ints
        All valid zip codes as integers. Zero-pad to length 5 to
        obtain the full zip code: `ZIP_string = "%05d" % zip`
    """

    # The file with all ZIP codes in the country in question. These
    # will be used to check the valitidity of user ZIP codes.
    zipDataFileName = ("/proj/net_scratch/data/MobilePhoneData/"
                       "GeoPostCode/GeoPC_ES/GeoPC_ES.csv")
    zips = set()
    with open(zipDataFileName, 'r') as f:
        for line in f:
            zips.add(int(line.split(';')[8][1:-1]))
            #zips.add(int(line.split("|")[0]))
    return zips

class UserData(np.recarray):
    """Represent user information.

    Note that differently from the original data file, missing,
    invalid and erroneous values are always marked with 0. This allows
    natural tests like `if gender` or `if usertype`.

    See the documentation of this file for more information on
    available user data fields.

    Usage:
    >>> from phone.userdata import UserData
    >>> ud = UserData('userData_mutual.txt')
    >>> userID = 1127
    >>> print ud[userID]['age'], ud[userID]['gender']
    52 1
    >>> print ud[userID]['ZIP'], ud[userID]['usertype']
    40297 0
    >>>
    >>> ud.save('userData_mutual.npy')
    >>> del ud
    >>> ud = UserData('userData_mutual.npy')
    """
    # Number of fields in a normal line.
    N_fields = 10

    # Define external constants.
    missing = unknown = 0
    prepaid, postpaid = 1, 2
    female, male = 1, 2
    max_age = 150

    # Define values used in the data.
    __invalid_ZIP_in_data = "99999"
    __invalid_age_in_data = 0
    __invalid_gender_in_data = 0
    __usertype_map = {'0': postpaid, '1': prepaid}
    __gender_map = {'0': unknown, '1': female, '2': male}

    def __initFromASCII(self, fileName):
        """Read user data from ASCII file.

        The columns of the file are
            [0] User ID (integer from 0 to N)
            [1] Age (restricted to [1, 150], 0 means unknown)
            [2] Gender (0 = unknown, 1 = female, 2 = male)
            [3] ZIP code (99999 is unknown)
            [4] User type (0 = postpaid, 1 = prepaid, 2 = unknown)
            [5] Contract ID
            [6] Activation date (YYYYMMDD)
            [7] Disconnection date (YYYYMMDD), 0 if not present
            [8] Latitude of most common position (0 = unknown)
            [9] Longitude of most common position (0 = unknown)
        See the documentation of this class for rules applied when
        interpreting missing or invalid values.
        """
        # Save the original file name.
        self.__originalASCIIFileName = fileName

        # Get a set of all valid ZIP codes.
        valid_zips = getValidZIPs()

        with open(fileName, 'rU') as f:
            for line in f:
                fields = line.split()
                userID = int(fields[0])
                age = min(UserData.max_age, max(int(fields[1]), 0))

                if ((len(fields) == 9 and "  " in line
                     and fields[3] in ("0","1"))
                    or len(fields[3]) != 5):
                    # ZIP code is missing in data or does not consist
                    # of exactly five characters => ZIP is corrupted.
                    ZIP_code = int(UserData.__invalid_ZIP_in_data)
                else:
                    # fields[3] has five characters. If converting to
                    # integer succeeds and that integer is included in
                    # the list of all valid integers, this is a valid
                    # ZIP code.
                    try:
                        ZIP_code = int(fields[3])
                    except ValueError:
                        ZIP_code = int(UserData.__invalid_ZIP_in_data)

                validZIP = False
                if ZIP_code in valid_zips:
                    validZIP = True


                # Because ZIP code (field 3) might be missing, we use
                # negative indexing for indices after it.
                usertype = UserData.__usertype_map.get(fields[-6],
                                                       UserData.unknown)

                gender = UserData.__gender_map[fields[2]]
                if (gender == UserData.male and 
                    usertype is not UserData.postpaid and 
                    (age == UserData.__invalid_age_in_data or 
                     fields[3] == UserData.__invalid_ZIP_in_data)):
                    gender = UserData.unknown

                contractID = int(fields[-5])
                activationDate = fields[-4]
                if activationDate == "0" or len(activationDate) != 8:
                    activationDate = ""
                disconnectDate = fields[-3]
                if disconnectDate == "0" or len(disconnectDate) != 8:
                    disconnectDate = ""
                latitude = float(fields[-2])
                longitude = float(fields[-1])

                # Add new record.
                self[userID]['age'] = age
                self[userID]['gender'] = gender
                self[userID]['ZIP'] = ZIP_code
                self[userID]['validZIP'] = validZIP
                self[userID]['usertype'] = usertype
                self[userID]['contractID'] = contractID
                self[userID]['activationDate'] = activationDate
                self[userID]['disconnectDate'] = disconnectDate
                self[userID]['latitude'] = latitude
                self[userID]['longitude'] = longitude
                      
    def __new__(cls, fileName, max_user_id=None):
        """Initialize user data from file.

        Parameters
        ----------
        fileName : str
            If the ending is '.npy', the file is assumed to be a saved
            UserData object and is loaded accordingly. Otherwise
            fileName must be a text file of user data.
        max_user_id : int
            The maximum user id. If None, the number of lines in
            `fileName` is used.
        """
        if fileName.split('.')[-1] == 'npy':
            # Read existing saved data and return the read object.
            return np.load(fileName)
        else:
            # Attempt to read text file. If the number of users is not
            # specified, the number of lines in fileName is used
            # instead.
            if max_user_id is None:
                max_user_id = int(os.popen("wc -l %s" % fileName).next().split()[0])
                print "Found %d lines." % max_user_id # DEBUG
                max_user_id
            # Construct the data object.
            dtype = [('age', np.uint8), ('gender', np.uint8), ('ZIP', np.uint32), 
                     ('validZIP', bool), ('usertype', np.uint8),
                     ('contractID', np.uint32), ('activationDate', '|S8'), 
                     ('disconnectDate', '|S8'), ('latitude', np.float32), 
                     ('longitude', np.float32)]
            obj = super(UserData, cls).__new__(cls, max_user_id+1, dtype=dtype)
            # Read in the data into the object.
            obj.__initFromASCII(fileName)
            return obj

    def originalFileName(self):
        """Returns the file name of the original data file."""
        return self.__originalASCIIFileName

    def save(self, fileName):
        """Save user data as numpy binary.

        File name must end with '.npy', because this is used to
        identify the file format when reading the file.
        """
        if fileName.split('.')[-1] != 'npy':
            raise ValueError("File name must end with '.npy'.")
        np.save(fileName, self)

    def isEmpty(self, userRecord):
        """Return True if data record of `userID` is empty.

        Because all users must have a non-zero contract id, this is
        the only field that needs to be checked.
        """
        return (userRecord['contractID'] == 0)

    def saveAsASCII(self, fileName):
        """Save user data as ASCII.

        Missing and non-numeral ZIP codes will be replaced by 99999.
        """
        if fileName.split('.')[-1] != 'txt':
            raise ValueError("File name must end with '.txt'.")
        formatStr = "%d %d %d %05d %d %d %s %s %s %s\n"
        genderMap = {UserData.unknown: 0, UserData.female: 1,
                     UserData.male: 2}
        typeMap = {UserData.unknown: 2, UserData.postpaid: 0,
                   UserData.prepaid: 1}
        with open(fileName, 'w') as f:
            for i, u in enumerate(self):
                if not self.isEmpty(u):
                    f.write(formatStr % (i, u['age'], genderMap[u['gender']],
                                         u['ZIP'], typeMap[u['usertype']],
                                         u['contractID'], u['activationDate'],
                                         (u['disconnectDate'] or "0"),
                                         str(u['latitude']), str(u['longitude'])))


if __name__ == '__main__':
    #from tests.test_userdata2 import *
    #unittest.main()
    pass
