# Copyright (C) 2016 by Per Unneberg
# Helper functions for filetypes
import os
try:
    from pytest_ngsfixtures import ROOT_DIR
except:
    ROOT_DIR = os.curdir

MAPPING = {
    'bai': "PUR.HG00731.tiny.sort.bai",
    'bam': "PUR.HG00731.tiny.sort.bam",
}
    
def filetype_mapping(ft, end="pe"):
    """Retrieve mapping for filetype to fixture file name.

    Params:
      ft (str): filetype
      end (str): sequencing mode

    Returns:
      str: fixture file names
    """
    if ft is None:
        return None
    try:
        fn = os.path.join(ROOT_DIR, "data", "applications", end, MAPPING[ft])
    except Exception as e:
        raise e
    return fn

