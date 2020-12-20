import os
from os import path
import json
import sys
from typing import Dict, List

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility import logger
    from Bio import SeqIO
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

mitoflex_dir = path.abspath(path.join(path.dirname(__file__), '..'))
profile_dir = path.abspath(path.join(mitoflex_dir, 'profile'))
genetic_code_profile = path.abspath(path.join(profile_dir, 'codes.json'))
hmm_dir = path.abspath(path.join(profile_dir, 'CDS_HMM'))
mtdb_dir = path.abspath(path.join(profile_dir, 'MT_database'))
required_cds_profile = path.abspath(path.join(hmm_dir, 'required_cds.json'))

if any(not path.isdir(x) for x in [mitoflex_dir, profile_dir, hmm_dir, mtdb_dir]) \
        or not path.isfile(genetic_code_profile):
    sys.exit(
        f"Profile structural check failed, is the installation of MitoFlex valid?")


def check_clade(clade_name: str) -> bool:
    return path.isfile(path.join(hmm_dir, f'{clade_name}.hmm')) or path.isfile(path.join(mtdb_dir, f'{clade_name}.fa'))


def delete_clade(clade_name: str):
    # Removing existed profiles
    if path.isfile(path.join(hmm_dir, f'{clade_name}.hmm')):
        logger.log(2, "Removing hmm profile.")
        os.remove(path.join(hmm_dir, f'{clade_name}.hmm'))
    if path.isfile(path.join(mtdb_dir, f'{clade_name}.fa')):
        logger.log(2, "Removing mt protein profile.")
        os.remove(path.join(mtdb_dir, f'{clade_name}.fa'))

    # Removing things in *.jsons.
    required_cds: dict = json.load(open(required_cds_profile))
    if clade_name in required_cds:
        logger.log(2, "Removing cds requirement in profile.")
        required_cds.pop(clade_name)
    json.dump(required_cds, open(required_cds_profile, 'w'))

    codes: Dict[str, int] = json.load(open(genetic_code_profile))
    if clade_name in codes:
        logger.log(2, "Removing genetic code profile.")
        codes.pop(clade_name)
    json.dump(codes, open(genetic_code_profile, 'w'))


def config_clade(clade_name: str, genetic_code: int, hmm_profile: str, mtdb_profile: str, cds_profile: str):
    # Dumping genetic code
    codes: Dict[str, int] = json.load(open(genetic_code_profile))
    codes[clade_name] = genetic_code
    json.dump(codes, open(genetic_code_profile))

    # Loading cds profile, and checking file integrity.
    to_insert: Dict[str, int] = json.load(open(cds_profile))
    