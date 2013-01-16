#!/usr/bin/env python2
from __future__ import print_function

import os
import re

def add_library_calls(libs, dirname, names):
    for name in names:
        if name.lower().endswith(".r"):
            path = os.path.join(dirname, name)
            for lib in re.findall("library\((.+?)\)", open(path).read()):
                if "=" in lib:
                    continue
                lib = lib.replace('"','')
                libs.add(lib)


libs = {"org.Hs.eg.db", "org.Mm.eg.db", "hom.Hs.inp.db"}
os.path.walk(".", add_library_calls, libs)
print("source(\"http://bioconductor.org/biocLite.R\")")
print(*['try({biocLite("'+lib+'")})' for lib in sorted(libs)],
      sep="\n")
