#!/usr/bin/env python2.7
import unittest
# append parent directory to import path
import env
# now we can import the lib module
from src import count_spike_in_seq

class TestSeqSubstring(unittest.TestCase):
    def test_long(self):
        self.assertEqual(count_spike_in_seq.get_ss_check_seq('1234567890ABCDEFGHIJKLMNO', 5, 10), 
                         '67890ABCDE')

    def test_short(self):
        self.assertEqual(count_spike_in_seq.get_ss_check_seq('1234567890ABCDEFGHIJKLMNO', 20, 10), 
                         'KLMNO')

    def test_none(self):
        self.assertEqual(count_spike_in_seq.get_ss_check_seq('1234567890ABCDEFGHIJKLMNO', 25, 10), 
                         '')


if __name__ == '__main__':
    unittest.main()