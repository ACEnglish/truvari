# Pull the current version of the wiki as a VERSION of the docs
VERSION=$1

set -e

if [ -e $1 ]
then
    echo 'Must provide version argument'
    exit
fi

git clone https://github.com/spiralgenetics/truvari.wiki.git
mv truvari.wiki $VERSION

echo 'Review, git add, and git push'
