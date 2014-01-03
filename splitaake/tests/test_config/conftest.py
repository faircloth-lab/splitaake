import os
import ConfigParser
import pytest

from splitaake import config


@pytest.fixture(scope="module")
def lq():
    return config.ListQueue()


# setup conf for config file
@pytest.fixture(scope="module")
def conf_file():
    conf = ConfigParser.ConfigParser()
    conf.read('splitaake-config-test.conf')
    return conf


@pytest.fixture(scope="module")
def parallelism(conf_file):
    return config.ConfParallelism(conf_file)


@pytest.fixture(scope="module")
def db(conf_file):
    return config.ConfDb(conf_file)


@pytest.fixture(scope="module")
def reads(conf_file):
    return config.ConfReads(conf_file)


@pytest.fixture(scope="module")
def quality(conf_file):
    return config.ConfQuality(conf_file)
