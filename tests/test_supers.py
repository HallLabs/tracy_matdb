# """Tests the supercell matrix selector that uses `supercell`.
# """
# import pytest
# from matdb.atoms import Atoms
# from matdb.transforms import _get_supers
# from matdb.utility import reporoot
# from os import path

# def test_supers():
#     at = Atoms(path.join(reporoot, "tests", "files", "POSCAR-AgPd-50"),
#                format="vasp")
#     sizes = [8, 11, 17]
#     result = _get_supers(at, sizes)

#     for s in sizes:
#         assert s in result
#     assert result[8].size == 8
#     assert result[11].size == 12
#     assert result[17].size == 16
