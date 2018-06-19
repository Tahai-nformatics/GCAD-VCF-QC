from unittest import TestCase
from gcad_vcf_qc import count_gt1


class TestCount_gt1(TestCase):

    dict1 = ([[('GT', (0, 0)), ('GQ', 54), ('DP', 7), ('HQ', (56, 60))],
              [('GT', (0, 0)), ('GQ', 48),
               ('DP', 4), ('HQ', (51, 51))],
              [('GT', (0, 0)), ('GQ', 61),
               ('DP', 2), ('HQ', (None,))],
              [('GT', (0, 0)), ('GQ', 48),
               ('DP', 1), ('HQ', (51, 51))],
              [('GT', (1, 0)), ('GQ', 48),
               ('DP', 8), ('HQ', (51, 51))],
              [('GT', (1, 1)), ('GQ', 43),
               ('DP', 5), ('HQ', (None, None))],
              [('GT', (0, 0)), ('GQ', 49),
               ('DP', 3), ('HQ', (58, 50))],
              [('GT', (0, 1)), ('GQ', 3),
               ('DP', 5), ('HQ', (65, 3))],
              [('GT', (0, 0)), ('GQ', 41),
               ('DP', 3), ('HQ', (None,))],
              [('GT', (1, 2)), ('GQ', 21),
               ('DP', 6), ('HQ', (23, 27))],
              [('GT', (2, 1)), ('GQ', 2),
               ('DP', 0), ('HQ', (18, 2))],
              [('GT', (2, 2)), ('GQ', 35),
               ('DP', 4), ('HQ', (None,))],
              [('GT', (0, 1)), ('GQ', 35), ('DP', 4)],
              [('GT', (0, 2)), ('GQ', 17), ('DP', 2)],
              [('GT', (1, 1)), ('GQ', 40), ('DP', 3)]])

    def test_count_gt1(self):

        dict2 = ([3, 1, 2, 0, 0])
        self.assertEqual(count_gt1(self.dict1), dict2)

    def test_demo(self):
        self.fail()
