"""
Module for analysing mobile phone network metadata.

NOTE: This file is used for mobile phone data set 1. For set 2, use
userdata2.py instead.

Usage:

>>> # Read user data either from a text file or from a previously
>>> # saved npy-file. Reading from npy-file is _much_ faster!
>>> from phone.userdata import UserData
>>> udata = UserData('userdata.txt')
>>> udata.save('userdata.npy')
>>> del udata
>>> udata = UserData('userdata.npy')

>>> # Take a look at the data of a single user.
>>> udata[1127]
(52, 1, 40297L, False, 2)
>>> # The previous line means that user 1127 is 52 year old female
>>> # (the second '1'), with valid ZIP 40297 (invalidZIP = False), and
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
>>> if not udata[1127]['invalidZIP']:
...     print "ZIP %05d is valid!" % udata[1127]['ZIP']
...
ZIP 40297 is valid!

>>> # Type UserData.unknown evaluates to False. Therefore it is
>>> # possible to go through valid users like this:
>>> N_valid_users = 0
>>> for user in udata:
...     if (user['gender'] and user['usertype'] and not user['invalidZIP']):
...         N_valid_users += 1
...
>>> print "Only %d valid users out of %d!" % (N_valid_users, len(udata))
Only 4162538 valid users out of 5343749!

"""

import os
import numpy as np

class UserData(np.recarray):
    """Represent user information.

    Note that differently from the original data file, missing,
    invalid and erroneous values are always marked with 0. This allows
    natural tests like `if gender` or `if usertype`. The following
    five fields exist for each user:

    age         The age of the user in years. Missing or invalid value
                is indicated with 0.

    gender      0 = Unknown, 1 = Female, 2 = Male. Gender will be
                marked as unknown for non-postpaid users if in the
                original data gender is male and either 1) ZIP is
                '99999' or 2) age is '0'.

    ZIP         The ZIP code of the user, missing value is 99999. Note
                that the ZIP code is saved as integer for performance
                reasons. Because ZIP codes always have exactly 5
                digits, you need to zero-pad the integer to get the
                original ZIP code, i.e. `ZIP_string = "%05d" % ZIP`

    invalidZIP  True if ZIP code is clearly invalid. This happens if
                1) The ZIP in the original data is '99999'.
                2) The ZIP in the original data contains letters 
                   (e.g. '34OO1').
                3) The ZIP in the original data has less than 5
                   digits.
                4) The first two digits of the ZIP are not between
                   '01' and '52', as required by the standard.
                Note that a ZIP code can still be non-existant even if
                invalidZIP is False; however, it is always invalid if
                invalidZIP is True.
                
    usertype    0 = unknown, 1 = prepaid, 2 = postpaid. If a line in
                the original file contains only 4 columns, usertype
                will be unknown.

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
    # Define external constants.
    missing = unknown = 0
    prepaid, postpaid = 1, 2
    female, male = 1, 2
    invalid_ZIP = 99999
    max_age = 125

    # Define values used in the data.
    __invalid_ZIP_in_data = "99999"
    __invalid_age_in_data = 0
    __invalid_gender_in_data = 0
    __usertype_map = {'0': postpaid, '1': prepaid}

    def __initFromASCII(self, fileName):
        """Read user data from ASCII file.

        The columns of the file are
            [0] User ID (integer from 0 to N)
            [1] Age (restricted to [1, 125], 0 means unknown)
            [2] Gender (0 = unknown, 1 = female, 2 = male)
            [3] ZIP code (99999 is unknown)
            [4] User type (0 = postpaid, 1 = prepaid, 2 = unknown)
            
        See the documentation of this class for rules applied when
        interpreting missing or invalid values.
        """
        with open(fileName, 'rU') as f:
            for line in f:
                fields = line.split()
                userID = int(fields[0])
                age = min(UserData.max_age, max(int(fields[1]), 0))

                if len(fields) < 5:
                    usertype = UserData.unknown
                else:
                    usertype = UserData.__usertype_map.get(fields[4],
                                                           UserData.unknown)

                gender = int(fields[2])
                if (gender == UserData.male and 
                    usertype is not UserData.postpaid and 
                    (age == UserData.__invalid_age_in_data or 
                     fields[3] == UserData.__invalid_ZIP_in_data)):
                    gender = UserData.unknown

                invalidZIP = True
                try:
                    ZIP_code = int(fields[3])
                    area_code = int(fields[3][:2])
                    if area_code > 0 and area_code < 53 and len(fields[3]) == 5:
                        invalidZIP = False
                except ValueError:
                    ZIP_code = UserData.invalid_ZIP
                

                # Add new record.
                self[userID]['age'] = age
                self[userID]['gender'] = gender
                self[userID]['ZIP'] = ZIP_code
                self[userID]['invalidZIP'] = invalidZIP
                self[userID]['usertype'] = usertype
                      
    def __new__(cls, fileName, N_users=None):
        """Initialize user data from file.

        Parameters
        ----------
        fileName : str
            If the ending is '.npy', the file is assumed to be a saved
            UserData object and is loaded accordingly. Otherwise
            fileName must be a text file of user data.
        N_users : int
            The number of users. If None, the number of lines in
            fileName is used.
        """
        if fileName.split('.')[-1] == 'npy':
            # Read existing saved data and return the read object.
            return np.load(fileName)
        else:
            # Attempt to read text file. If the number of users is not
            # specified, the number of lines in fileName is used
            # instead.
            if N_users is None:
                N_users = int(os.popen("wc -l %s" % fileName).next().split()[0])
                #print "Found %d lines." % N_users # DEBUG
            # Construct the object.
            dtype = [('age', np.uint8), ('gender', np.uint8), ('ZIP', np.uint32), 
                     ('invalidZIP', bool), ('usertype', np.uint8)]
            obj = super(UserData, cls).__new__(cls, N_users, dtype=dtype)
            # Read in the data into the object.
            obj.__initFromASCII(fileName)
            return obj

    def save(self, fileName):
        """Save user data as numpy binary.

        File name must end with '.npy', because this is used to
        identify the file format when reading the file.
        """
        if fileName.split('.')[-1] != 'npy':
            raise ValueError("File name must end with '.npy'.")
        np.save(fileName, self)

if __name__ == '__main__':
    from tests.test_userdata import *
    unittest.main()
