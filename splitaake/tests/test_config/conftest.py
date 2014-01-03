import os
import ConfigParser
import pytest

from splitaake import config


@pytest.fixture(scope="module")
def lq():
    return config.ListQueue()


@pytest.fixture(scope="module")
def cwd():
    return os.path.dirname(os.path.realpath(__file__))


# setup conf for config file
@pytest.fixture(scope="module")
def conffile():
    conf = ConfigParser.ConfigParser()
    conf.read('splitaake-config-test.conf')
    return conf


@pytest.fixture(scope="module")
def tags(conffile):
    return config.Tags(conffile)


@pytest.fixture(scope="module")
def test_tags():
    p = config.Tags(False, 0, 'forward')
    combos = (
        ("p1,n1", "combo1"),
        ("p1,n2", "combo2"),
        ("n1,p1", "combo3")
    )
    p._get_combinations(combos)
    tags = (
        ('p1', 'TTTT'),
        ('n1', 'GGGG'),
        ('n2', 'AAAA')
    )
    p._get_sequences(tags)
    return p
