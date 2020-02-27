import unittest
import os
if __name__ == '__main__':
    import phone.userdata as userdata
else:
    import userdata

class Test_UserData(unittest.TestCase):
    def setUp(self):
        self.testFileName = "test_userdata.%d.txt" % os.getpid()
        
        # Test data. The first vector will be written into the file,
        # the second one is the expected content of UserData.
        unk = userdata.UserData.unknown
        post = userdata.UserData.postpaid
        pre = userdata.UserData.prepaid
        self.ages = ([200, 0, 78, 55, 13, 100],
                     [userdata.UserData.max_age, 0, 78, 55, 13, 100])
        self.genders = ([1, 0, 1, 2, 1, 2],
                        [1, unk, 1, 2, 1, unk])
        self.ZIPs = (['01000', '34852', '56o67', '53110', '00100', '99999'],
                     ['01000', '34852', '99999', '53110', '00100', '99999'])
        self.invalidZIPs = (False, False, True, True, True, True)
        self.userTypes = ([0, 0, 1, 1, 2, ''],
                          [post, post, pre, pre, unk, unk])
        with open(self.testFileName, 'w') as f:
            userID = 0
            for i, ud in enumerate(zip(self.ages[0], self.genders[0],
                                       self.ZIPs[0], self.userTypes[0])):
                newLine = " ".join(map(str, [i] + list(ud)))
                f.write(newLine + "\n")

    def tearDown(self):
        os.remove(self.testFileName)

    def test_construction(self):
        """Create user data from text file."""
        ud = userdata.UserData(self.testFileName)

        # Make sure the data is good.
        for i, u in enumerate(ud):
            self.assertEqual(u.age, self.ages[1][i])
            self.assertEqual(u.gender, self.genders[1][i])
            self.assertEqual("%05d" % u.ZIP, self.ZIPs[1][i])
            self.assertEqual(u.invalidZIP, self.invalidZIPs[i])
            self.assertEqual(u.userType, self.userTypes[1][i])
        
    def test_save(self):
        """Save data as binary and restore."""
        ud_a = userdata.UserData(self.testFileName)
        npyFileName = self.testFileName[:-3] + 'npy'
        ud_a.save(npyFileName)
        ud_b = userdata.UserData(npyFileName)

        # Make sure ud_a and ud_b are identical.
        for ua, ub in zip(ud_a, ud_b):
            self.assertEqual(ua, ub)

        # Remove the binary data file.
        os.remove(npyFileName)

    def test_single(self):
        """Get just one data field."""
        ud = userdata.UserData(self.testFileName)
        ud[0].age

if __name__ == '__main__':
    #suite = unittest.TestSuite()
    #suite.addTest(Test_UserData("test_single"))
    #unittest.TextTestRunner().run(suite)
    unittest.main()
