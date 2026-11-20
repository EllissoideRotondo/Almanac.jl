"""Tests of signatures.py on small Fortran sources written for the tests."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import signatures  # noqa: E402

SUB = """\
C$Procedure SPKX ( A test routine )

      SUBROUTINE SPKX ( TARG, STARG,
     .                  X, A, B )

      IMPLICIT NONE
      INCLUDE 'test.inc'

      CHARACTER*(*)         TARG
      DOUBLE PRECISION      STARG    ( 6 )
      INTEGER               X
      DOUBLE PRECISION      A        ( 3, NDIM )
      INTEGER               B        ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     TARG       I   Target body name.
C     STARG      O   State of target, a description that goes on
C                    to a second line.
C     X         I-O  Something updated.
C     A,
C     B          I   Two inputs on two lines.
C
C$ Detailed_Input
C
C     TARG     is the name.
C
C-&
      INTEGER               LOCAL ( 2 )
      X = X + 1
      RETURN
      END
"""

UMBRELLA = """\
C$Procedure UMB ( Umbrella )

      SUBROUTINE UMB ( A, B, C )

      IMPLICIT NONE
      INTEGER               A
      DOUBLE PRECISION      B ( 3 )
      LOGICAL               C

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     A          I   Entry dependent.
C     B          I   Entry dependent.
C     C          O   Entry dependent.
C
C-&
      RETURN

C$Procedure ENT1 ( An entry point )

      ENTRY ENT1 ( A, B )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     A          I   An input.
C     B          O   An output.
C
C-&
      RETURN
      END
"""

FUNCTION = """\
C$Procedure VNX ( A function )

      DOUBLE PRECISION FUNCTION VNX ( V, WIN )

      IMPLICIT NONE
      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )
      DOUBLE PRECISION      V ( 3 )
      DOUBLE PRECISION      WIN ( LBCELL:* )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     V          I   A vector.
C     WIN       I/O  A window.
C
C     The function returns a number.
C
C-&
      VNX = V(1)
      RETURN
      END
"""

PROCS = """\
      SUBROUTINE WNX ( RESULT, BAIL )

      IMPLICIT NONE
      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )
      DOUBLE PRECISION      RESULT ( LBCELL:* )
      LOGICAL               BAIL
      EXTERNAL              BAIL

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     RESULT     O   A window.
C     BAIL       I   A function.
C
C-&
      RETURN
      END
"""

INCLUDES = {"test.inc": "      INTEGER NDIM\n      PARAMETER ( NDIM = 4 )\n"}


class SignatureTests(unittest.TestCase):
    def test_subroutine(self):
        t = signatures.parse_file(SUB, INCLUDES)["spkx"]
        self.assertEqual(t["kind"], "subroutine")
        self.assertEqual(t["ftype"], "")
        self.assertEqual(t["args"], [
            ("targ", "c", "", "in"),
            ("starg", "d", "6", "out"),
            ("x", "i", "", "inout"),
            ("a", "d", "3,4", "in"),
            ("b", "i", "*", "in"),
        ])

    def test_entry(self):
        t = signatures.parse_file(UMBRELLA, {})
        self.assertEqual(t["umb"]["args"], [("a", "i", "", "in"), ("b", "d", "3", "in"), ("c", "l", "", "out")])
        self.assertEqual(t["ent1"]["kind"], "entry")
        self.assertEqual(t["ent1"]["args"], [("a", "i", "", "in"), ("b", "d", "3", "out")])

    def test_function(self):
        t = signatures.parse_file(FUNCTION, {})["vnx"]
        self.assertEqual(t["kind"], "function")
        self.assertEqual(t["ftype"], "d")

    def test_cell_dims(self):
        t = signatures.parse_file(FUNCTION, {})["vnx"]
        self.assertEqual(t["args"][1], ("win", "d", "LBCELL:*", "inout"))

    def test_output_cell_is_read(self):
        t = signatures.parse_file(PROCS, {})["wnx"]
        self.assertEqual(t["args"][0], ("result", "d", "LBCELL:*", "inout"))

    def test_typed_procedure(self):
        t = signatures.parse_file(PROCS, {})["wnx"]
        self.assertEqual(t["args"][1], ("bail", "p", "", "in"))

    def test_round_trip(self):
        t = signatures.parse_file(SUB, INCLUDES)
        path = os.path.join(os.environ.get("TMPDIR", "/tmp"), "sig_test.tsv")
        signatures.write(t, path)
        self.assertEqual(signatures.load(path), t)


if __name__ == "__main__":
    unittest.main()
