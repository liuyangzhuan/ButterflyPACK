cat $1.tar.gz_* > $1.tar.gz
rm -rf $1.tar.gz_*
tar -xvf $1.tar.gz
rm -rf $1.tar.gz

