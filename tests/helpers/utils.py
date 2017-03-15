# Copyright (C) 2016 by Per Unneberg
# Helper function to make output executable
import os

def _make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

def save_command(fn, args):
    with open(fn, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write("PATH={}\n".format(os.environ["PATH"]))
        fh.write("args=$*\n")
        fh.write(" ".join(args) + " ${args}\n")
    _make_executable(fn)

