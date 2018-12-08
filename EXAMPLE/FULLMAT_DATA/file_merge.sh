cat $1.tar.gz_* > $1.tar.gz
rm $1.tar.gz_*
tar -xvf $1.tar.gz
rm $1.tar.gz

