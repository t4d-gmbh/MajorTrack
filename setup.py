from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
        name='MajorTrack',
        # version='0.1',
        # ####################################################################
        # ####################################################################
        use_scm_version={'write_to': 'docs/version.txt'},
        setup_requires=['setuptools_scm'],
        # ####################################################################
        # ####################################################################
        description='Evolutionary clustering method with fission-fusion '
        'detection.',
        long_description=readme(),
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'License :: OSI Approved :: GNU General Public License 3',
            'Programming Language :: Python :: 3.6',
            'Topic :: Evolutionary Clustering :: Machine Learning :: '
            'Dynamic Community',
        ],
        keywords='clustering community classification temporal',
        url='https://github.com/j-i-l/majortrack',
        author='Jonas I. Liechti [aut, cre]',
        author_email='jonas.i.liechti@gmail.com',
        license='GPL-3',
        packages=['majortrack'],
        install_requires=[
          'matplotlib', 'pyalluv==0.1', 'colorseq==0.1'
        ],
        dependency_links=[
            'git+https://github.com/j-i-l/pyalluv.git@v0.1#egg=pyalluv-0.1',
            'git+https://github.com/j-i-l/'
            'colorsequence.git@v0.1#egg=colorseq-0.1'
            ],
        test_suite='nose.collector',
        tests_require=['nose', 'nose-cover3'],
        # ToDo:
        # entry_points={
        #   'console_scripts': ['draw-alluvial=pyalluv.command_line:main'],
        # },
        include_package_data=True,
        zip_safe=False
        )
