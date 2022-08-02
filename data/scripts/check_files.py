#!/usr/bin/env python

_version__ = "0.1"
__author__ = "Michi Jackson"
__contact__ = "noe.pozzan@stud.unibas.ch"
__doc__ = "Parse a parameter file and change some parameters."

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import os
from pathlib import Path
import sys

def main():
    """ Main function """

	number = int(sys.argv[1])

	for i in range(2, number+2):
		file = Path(sys.argv[i])
		if file.exists():
			print("File {} is present".format(file))
		else:
			raise FileNotFoundError("No such file or directory: {}".format(file))



# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)


