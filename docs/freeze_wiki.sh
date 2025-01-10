# Pull the current version of the wiki 

git clone https://github.com/spiralgenetics/truvari.wiki.git
rm -rf truvari.wiki/.git
mv truvari.wiki/* ./
rm -rf truvari.wiki

echo 'Review, git add, and git push'
