from setuptools import setup, find_packages

setup(
    name='MemXTerminator',
    version='1.0.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    author='Zhen Huang',
    author_email='zhen.victor.huang@gmail.com',
    description='A software for membrane analysis and subtraction in cryo-EM',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    license='GPLv3',
    url='',
    keywords=['cryo-EM', 'cryo-ET', 'membrane subtraction'],
    python_requires='>=3.9',
    install_requires=[
        'cryosparc-tools>=4.3.1',
        'deap>=1.4.1',
        'matplotlib>=3.7.2',
        'mrcfile>=1.4.3',
        'multiprocess>=0.70.15',
        'numpy>=1.25.2',
        'pandas>=2.1.0',
        'PyQt5>=5.15.9',
        'scikit-image>=0.19.3',
        'scipy>=1.11.1',
        'starfile>=0.4.12',
        'cupy-cuda12x>=12.2.0'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': [
            'MemXTerminator=memxterminator.cli.main_cli:main',
        ]
    }
)
