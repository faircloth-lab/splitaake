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
def sites(conffile):
    return config.Sites(conffile)


@pytest.fixture(scope="module")
def test_tags():
    p = config.Tags(False, 0, 'forward', False)
    combos = (
        ("p1,n1", "combo1"),
        ("p1,n2", "combo2"),
        ("n1,p1", "combo3")
    )
    p._get_combinations(combos)
    tags = (
        ('p1', 'TTTTtt'),
        ('n1', 'GGGGgg'),
        ('n2', 'AAAAaa')
    )
    p._get_sequences(tags)
    return p


@pytest.fixture(scope="module")
def test_sites():
    return config.Sites(False, 0, 'TGCA', 'TA', auto=True)
