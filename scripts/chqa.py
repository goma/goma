#!/usr/bin/python
# Usage:
#      chqa.py <filename>
#
# Simple python script that retrieves and prints to stdout
# the qa records of exodus file <filename>
#
# Note: This requires the exodus python module to be in a
# directory listed in $PYTHONPATH and a dynamic exodus
# library (exodus.so) to be named in exodus.py.

import exodus
import sys

ef = exodus.exodus(sys.argv[1], 'r')
print(ef.get_qa_records())
ef.close()
