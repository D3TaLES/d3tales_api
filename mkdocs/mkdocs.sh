# using sphinx

version=0.2
author='Duke, Rebekah'
projectname='D<sup>3</sup>TaLES API'

whereisit="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $whereisit
cd ../
rm -rf docs doc
export DB_INFO_FILE=$PWD/db_info_ex.json

sphinx-apidoc -H "$projectname" -A "$author" -V $version --full -o ./doc ./d3tales_api tests
cp $whereisit/conf.py ./doc/
cp $whereisit/index.rst ./doc/
cp $whereisit/*.md ./doc/
cp $whereisit/*.png ./doc/

cd doc
make clean
make html

mv _build/html ../docs
cd ../
rm -rf doc
touch docs/.nojekyll
cd $whereisit


