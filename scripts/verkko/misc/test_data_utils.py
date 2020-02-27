import unittest
import numpy as np
import data_utils as du
import os

class Test_DataUtils(unittest.TestCase):
    def setUp(self):
        # Write input files for testing file reading.
        self.filenames = []
        self.columns = {}
        self.header_lines = {}

        # Write correct data
        fn = 'test_data_utils.data.%d.%d' % (len(self.filenames), os.getpid())
        self.filenames.append(fn)
        self.header_lines[fn] = 2
        self.columns[fn] = [[0, 1, 2, 3],
                            [1, 2, 4, 8],
                            [0.5, 0.5, 0.5, 0.5],
                            [0.8, 1.6, 2.4, 3.2]]

        # This data will miss one element in third row.
        fn = 'test_data_utils.data.%d.%d' % (len(self.filenames), os.getpid())
        self.filenames.append(fn)
        self.header_lines[fn] = 0
        self.columns[fn] = [[0, 1, 2, 3],
                            [1, 2, 4, 8],
                            [0.5, 0.5, 0.5, 0.5],
                            [0.8, 1.6, " ", 3.2]]

        # This data will have one float on the first row.
        fn = 'test_data_utils.data.%d.%d' % (len(self.filenames), os.getpid())
        self.filenames.append(fn)
        self.header_lines[fn] = 1
        self.columns[fn] = [[0.5, 1, 2, 3],
                            [1, 2, 4, 8],
                            [0.5, 0.5, 0.5, 0.5],
                            [0.8, 1.6, 2.4, 3.2]]

        # Write data files
        for fn in self.filenames:
            with open(fn, 'w') as f:
                for hl_num in range(self.header_lines[fn]):
                    f.write("Header_line %d\n" % (hl_num,))
                for line_num in range(len(self.columns[fn][0])):
                    fields = [self.columns[fn][i][line_num] 
                              for i in range(len(self.columns[fn]))]
                    f.write(" ".join(map(str, fields)) + "\n")

    def tearDown(self):
        # Remove the data files.
        for fn in self.filenames:
            os.remove(fn)

    def test_read_columns(self):
        """Test basic column reading."""
        fn = self.filenames[0]
        
        # Get the first two columns. Both columns are of type int.
        col0, col1 = du.read_columns(fn, self.header_lines[fn], (int, int))
        self.assertEqual(col0, self.columns[fn][0])
        self.assertEqual(col1, self.columns[fn][1])

        # Read columns 1 and 3 and convert them to type
        # int and float, respectively. Skip the two header rows.
        col1, col3 = du.read_columns(fn, self.header_lines[fn], (int, float), (1,3))
        self.assertEqual(col1, self.columns[fn][1])
        self.assertEqual(col3, self.columns[fn][3])

        # Read only the second column. Note the use of tuples!
        col1, = du.read_columns(fn, self.header_lines[fn], (int,), (1,))
        self.assertEqual(col1, self.columns[fn][1])

        # A missing column should raise IndexError.
        fn = self.filenames[1]
        self.assertRaises(IndexError, du.read_columns,
                          fn, self.header_lines[fn], (int,float), (1,3))

        # Finding a float when int is expected should raise
        # ValueError.
        fn = self.filenames[2]
        self.assertRaises(ValueError, du.read_columns,
                          fn, self.header_lines[fn], (int,int), (0,1))

    def test_read_columns_fun(self):
        """Test column reading with conversion function."""
        fn = self.filenames[0]
        
        # Get the first two columns. Both columns are of type int.
        def col_fun(fields):
            return int(fields[0]), int(fields[1])
        col0, col1 = du.read_columns(fn, self.header_lines[fn], col_fun)
        self.assertEqual(col0, self.columns[fn][0])
        self.assertEqual(col1, self.columns[fn][1])

        # Read columns 1 and 3 and convert them to type
        # int and float, respectively. Skip the two header rows.
        def col_fun(fields):
            return int(fields[1]), float(fields[3])
        col1, col3 = du.read_columns(fn, self.header_lines[fn], col_fun)
        self.assertEqual(col1, self.columns[fn][1])
        self.assertEqual(col3, self.columns[fn][3])

        # Read only the second column. Note the use of tuples!
        def col_fun(fields):
            return int(fields[1]),
        col1, = du.read_columns(fn, self.header_lines[fn], col_fun)
        self.assertEqual(col1, self.columns[fn][1])

    def test_gen_columns(self):
        """Test basic column generation."""
        fn = self.filenames[0]

        # Get the first two columns. Both columns are of type int.
        data_gen = du.gen_columns(fn, self.header_lines[fn], (int, int))
        for i, (c0, c1) in enumerate(data_gen):
            self.assertEqual(c0, self.columns[fn][0][i])
            self.assertEqual(c1, self.columns[fn][1][i])

        # Read columns 1 and 3 and convert them to type
        # int and float, respectively. Skip the two header rows.
        data_gen  = du.gen_columns(fn, self.header_lines[fn], (int, float), (1,3))
        for i, (c1, c3) in enumerate(data_gen):
            self.assertEqual(c1, self.columns[fn][1][i])
            self.assertEqual(c3, self.columns[fn][3][i])

        # Read only the second column. Note the use of tuples!
        data_gen = du.gen_columns(fn, self.header_lines[fn], (int,), (1,))
        for i, (c1,) in enumerate(data_gen):
            self.assertEqual(c1, self.columns[fn][1][i])

        # Finding a float when int is expected should raise
        # ValueError.
        fn = self.filenames[2]
        data_gen = du.gen_columns(fn, self.header_lines[fn], (int,int), (0,1))
        self.assertRaises(ValueError, lambda y: [x for x in y], data_gen)

        # A missing column should raise IndexError.
        fn = self.filenames[1]
        data_gen = du.gen_columns(fn, self.header_lines[fn], (int,float), (1,3))
        self.assertRaises(IndexError, lambda y: [x for x in y], data_gen)
        
    def test_gen_columns_fun(self):
        """Test column generation with conversion function."""
        fn = self.filenames[0]

        # Get the first two columns. Both columns are of type int.
        def col_fun(fields):
            return int(fields[0]), int(fields[1])
        data_gen = du.gen_columns(fn, self.header_lines[fn], col_fun)
        for i, (c0, c1) in enumerate(data_gen):
            self.assertEqual(c0, self.columns[fn][0][i])
            self.assertEqual(c1, self.columns[fn][1][i])

        # Read columns 1 and 3 and convert them to type
        # int and float, respectively. Skip the two header rows.
        def col_fun(fields):
            return int(fields[1]), float(fields[3])
        data_gen  = du.gen_columns(fn, self.header_lines[fn], col_fun)
        for i, (c1, c3) in enumerate(data_gen):
            self.assertEqual(c1, self.columns[fn][1][i])
            self.assertEqual(c3, self.columns[fn][3][i])

        # Read only the second column. Note the use of tuples!
        def col_fun(fields):
            return int(fields[1]),
        data_gen = du.gen_columns(fn, self.header_lines[fn], col_fun)
        for i, (c1,) in enumerate(data_gen):
            self.assertEqual(c1, self.columns[fn][1][i])

    def test_cumulative_dist(self):
        # Function for converting integers to strings with given
        # accuracy.
        to_str = lambda x: "%.10f" % x

        # TEST DESCENDING: P(X >= x)
        correct_data = [0.1,   1,   4,   8,   9]
        correct_prob = [  1, 0.8, 0.4, 0.2, 0.1]
        correct_prob = map(to_str, correct_prob)

        # Test with data only
        data = [1, 0.1, 9, 4, 1, 1, 8, 1, 4, 1.0/10]
        cum_data, cum_prob = du.cumulative_dist(data)
        self.assertEqual(cum_data, correct_data)
        self.assertEqual(map(to_str, cum_prob), correct_prob)

        # Test with count
        data =  [1, 0.1, 9, 4, 8]
        count = [4,   2, 1, 2, 1]
        cum_data, cum_prob = du.cumulative_dist(data, count)
        self.assertEqual(cum_data, correct_data)
        self.assertEqual(map(to_str, cum_prob), correct_prob)
        
        # TEST ASCENTING: P(X <= x)
        correct_data = [0.1,   1,   4,   8,   9]
        correct_prob = [0.2, 0.6, 0.8, 0.9, 1.0]
        correct_prob = map(to_str, correct_prob)

        # Test with data only
        data = [1, 0.1, 9, 4, 1, 1, 8, 1, 4, 1.0/10]
        cum_data, cum_prob = du.cumulative_dist(data, format='ascenting')
        self.assertEqual(cum_data, correct_data)
        self.assertEqual(map(to_str, cum_prob), correct_prob)

        # Test with count
        data =  [1, 0.1, 9, 4, 8]
        count = [4,   2, 1, 2, 1]
        cum_data, cum_prob = du.cumulative_dist(data, count, format='ascenting')
        self.assertEqual(cum_data, correct_data)
        self.assertEqual(map(to_str, cum_prob), correct_prob)

        # TEST EXCEPTIONS
        bad_format = 'descenting'
        self.assertRaises(ValueError, du.cumulative_dist, data, format=bad_format)
        bad_count = [4,   1, 1, 2, 1, 7]
        self.assertRaises(ValueError, du.cumulative_dist, data, bad_count)


if __name__ == '__main__':
    if False:
        suite = unittest.TestSuite()
        suite.addTest(Test_DataUtils("test_gen_columns"))
        unittest.TextTestRunner().run(suite)
    else:
        unittest.main()
