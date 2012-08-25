import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup
from setuptools import find_packages
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension

if __name__ == '__main__':
    setup(
        name='splitaake',
        version="1.0",
        description="Demultiplex hierarchically sequence-tagged massively"
            +"parallel sequencing reads",
        author="Brant Faircloth",
        author_email="brant.faircloth+demuxipy@gmail.com ",
        url="http://github.com/faircloth-lab/demuxipy",
        license="http://www.opensource.org/licenses/BSD-3-Clause",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        requires=['numpy(>=1.3)','seqtools(>=0.5)'],
        long_description=open('README.rst').read(),
        scripts=['bin/demuxi.py',
                ],
        ext_modules=[
                Extension(
                    'demuxipy/cpairwise2',
                        [
                            'demuxipy/cpairwise2module.c'
                        ]
                    ),
                ],
        packages=find_packages(),
        include_package_data = True,
        )
