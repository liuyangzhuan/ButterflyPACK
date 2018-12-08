tar -cvzf $1.tar.gz $1
rm $1
split -b 40M $1.tar.gz $1.tar.gz_
rm $1.tar.gz
