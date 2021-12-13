"""
Example unit tests for drp_tools package
"""
import unittest
import desc.drp_tools

class drp_toolsTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.drp_tools.drp_tools(self.message)
        self.assertEqual(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.drp_tools.drp_tools)
        foo = desc.drp_tools.drp_tools(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
