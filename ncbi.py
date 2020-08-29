#!/usr/bin/env python3

"""
ncbi.py
========

Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

This file is part of MitoFlex.

MitoFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoFlex.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
from os import path


import sys
if sys.version_info[0] < 3:
    sys.exit('Python 3 must be installed in current environment! Please check if your environment setup (like conda environment) is deactivated or wrong!')

try:
    from ete3 import NCBITaxa
    try:
        import ast
        import inspect
        print("Patching NCBITaxa's base methods. For reason, see #2.\n")
        code_to_patch = """db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""
        patched_code = """db.execute("INSERT OR REPLACE INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""

        ncbiquery = sys.modules[NCBITaxa.__module__]
        lines_code = [x.replace(code_to_patch, patched_code)
                      for x in inspect.getsourcelines(ncbiquery.upload_data)[0]]
        lines_code.insert(1, "    print('\\nIf this message shown, then the patch is successful!')\n")
        lines_code.insert(1, "    import os, sqlite3, sys\n")
        lines_code.insert(1, "    DB_VERSION = 2\n")
        lines_code = "".join(lines_code)

        ast_tree = ast.parse(lines_code)
        patched_function = compile(ast_tree, "<string>", mode="exec")
        mod_dummy = {}
        exec(patched_function, mod_dummy)
        ncbiquery.upload_data = mod_dummy["upload_data"]
    except Exception:
        print("Patching failed, current taxonomy data downloaded from FTP may be failed to parse with ETE3!")
    finally:
        print("Patch finished.")

    import shutil
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

tot, used, free = shutil.disk_usage('/')
print("Filesystem status:",
      f"Total: {tot//(2**30):.2f} GB",
      f"Free: {free//(2**30):.2f} GB",
      "",
      "If the free disk space is too low (<1G), database updating can be failed!",
      sep='\n')


dump_file = path.join(path.dirname(__file__), 'taxdump.tar.gz')
dump_file = path.abspath(dump_file)
dump_dir = path.dirname(dump_file)
dump_file_old = path.join(dump_dir, 'old.taxdump.tar.gz')

use_old = ""
while use_old.upper() != "Y" and use_old.upper() != "N":
    use_old = input("Do you want to download database from NCBI, or use already-downloaded old taxonomy data? [y/n]:").strip()

if use_old.upper() == "Y":
    if os.path.isfile(dump_file):
        os.rename(dump_file, dump_file_old)

    try:
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        if os.path.isfile(dump_file_old):
            os.remove(dump_file_old)

        print("Testing database...")
        result = ncbi.get_name_translator(["Platyhelminthes"])  # May change to whatever needed
        if not result:
            print("Cannot retrieve taxid of Platyhelminthes, this may indicate a failure of database updating!")
        else:
            print("Successfully updated NCBI taxonomy database.")

    except Exception:
        print("Errors occured when fetching data from NCBI database, falling back to the last fetched database.")
        if path.isfile(dump_file_old):
            os.rename(dump_file_old, dump_file)
            ncbi = NCBITaxa(taxdump_file=os.path.abspath(dump_file))
        else:
            print("A taxdump file is not found under installation directory, cannot build NCBI taxanomy database.")
            print("Please manually download it from http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz , and move it to the installation directory.")
else:
    print("Using old taxonomy database, some new taxanomy entires may be missing.")
    ncbi = NCBITaxa(taxdump_file=os.path.abspath(dump_file))
