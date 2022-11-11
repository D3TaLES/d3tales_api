import setuptools
from d3tales_api import __version__, __author__,__credits__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='d3tales_api',
    version=__version__,
    author=__author__,
    author_email='d3tales@gmail.com',
    description='API for parsing, database connection, and property calculations for D3taLES project ',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/D3TaLES/d3tales_api',
    project_urls={
        "Bug Tracker": "https://github.com/D3TaLES/d3tales_api/issues"
    },
    license=__credits__,
    packages=['d3tales_api'],
    install_requires=['requests'],
)

