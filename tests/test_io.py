"""Tests the I/O functions for template hierarchy parsing in matdb.
"""
import pytest
from matdb.io import read
from matdb.utility import reporoot

def test_schema():
    """Tests the template read for the whole schema to make sure that contexts
    work as expected.
    """
    varfull = read(reporoot, "tests/io/schema")
    #The full dictionary is large and unwieldly and gets tested below in
    #separate pieces. For now, we just test that the context hookups worked
    #correctly. Edges uses a relative path specifier with `../` syntax, so it is
    #most likely to mess up.
    model = [
        {'targets': 'B',
         'name': 'AtoB',
         'doc': 'Connects A to B.',
         'sources': 'A',
         'properties': [
             {'dtype': 'int',
              'doc': 'some integer.',
              'example': 3,
              'name': 'value'}]
        }]
    assert varfull["edges"] == model

def test_raises_errors():
    """Makes sure errors are raised where appropriate.
    """
    with pytest.raises(ValueError):
        template = read(reporoot, "tests/dummy")

def test_corner_cases():
    """Tests the corner cases in the template reader.
    """
    varcorner = read(reporoot, "tests/io/corner")
    model = {
        'first': {'a': 0, 'b': 1},
        'second': [{'a': 0, 'b': 1}]
    }
    assert varcorner == model
        
def test_read_types():
    """Tests the recursive read of a directory of edge template files with
    sub-directories, relative paths, etc.
    """
    varkind = read(reporoot, "tests/io/types/varkind")
    model = {
        'name': 'VarKind',
        'properties': [
            {'name': 'tolerance',
             'example': 0.2,
             'dtype': 'float',
             'doc': 'heard of epsilon.'}
        ],
        'doc': 'Simple custom type for unit testing.'
    }
    assert varkind == model
        
def test_read_edges():
    """Tests the recursive read of a directory of edge template files with
    sub-directories, relative paths, etc.
    """
    AtoB = read(reporoot, "tests/io/edges/atob")
    model = {
        'targets': 'B',
        'properties': [
            {'dtype': 'int',
             'doc': 'some integer.',
             'name': 'value',
             'example': 3}],
        'doc': 'Connects A to B.',
        'name': 'AtoB',
        'sources': 'A'
    }
    assert AtoB == model
        
def test_read_verts():
    """Tests the recursive read of a directory of vertex template files with
    sub-directories, relative paths, etc.
    """
    A = read(reporoot, "tests/io/verts/a")
    modelA = {
        'name': 'A',
        'properties': [
            {'name': 'kind',
             'dtype': 'VarKind',
             'example': 'kind(1.8)',
             'keytype': 'float',
             'doc': 'testing how kind the vertex is.'
            }],
        'doc': 'Simple vertex for unit tests.'
    }
    assert A == modelA

    B = read(reporoot, "tests/io/verts/b")
    modelB = {
        'vtype': 'logic',
        'doc': 'Simple schema for unit tests.  ',
        'properties': [
            {'example': 3,
             'doc': 'some integer.',
             'dtype': 'int',
             'name': 'value'},
            {'example': "'cc'",
             'doc': 'b.text',
             'dtype': 'str',
             'name': 'text'}
        ],
        'name': 'B'}
    assert B == modelB
