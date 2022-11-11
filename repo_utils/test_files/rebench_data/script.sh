rm -rf test1/phab/
export PATH=`pwd`/fake_mafft/:$PATH
truvari rebench -u test1/ -f chr20.fa
