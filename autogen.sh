#!/bin/bash

#remove the old configure files.
rm -rf autom4te.cache configure config.* aclocal.* install-sh libtool ltmain.sh missing stamp-h1
rm -rf configure Makefile Makefile.in

aclocal
autoheader

# On MAC OS X, GNU libtoolize is named 'glibtoolize':
if [ `(uname -s) 2>/dev/null` == 'Darwin' ]
then
 echo "MAC glibtoolize"
 libtoolize --automake --copy
else
 echo "GNU libtoolize"
 libtoolize --automake --copy
fi 

automake --foreign --add-missing --copy
autoconf

