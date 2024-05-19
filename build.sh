# Prerequisites for slimmc:
#
# nim programming lang compiler 1.6.20
# https://nim-lang.org/install_unix.html
# https://nim-lang.org/install_windows.html
#
# gcc 10.1.*
# binutils
# git




GTAG=`git describe --abbrev=4 --dirty --always --tags`
GHASH=`git rev-parse HEAD`
BUILD=`date`
SYS=`uname`
NIM=`nim --version | head -n 1`
GCC=`gcc --version | head -n 1`

mkdir bin

echo '>>> Building slimmc'
PRG='slimmc'
EXTRA=''
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\nvar gcc=\"$GCC\"\n" > ./slimmc/version.nim
nim c -d:release -o:./bin/slimmc ./slimmc/slimmc.nim
rm ./slimmc/version.nim

echo ">>> Building slimmc-turbo"
PRG='slimmc-turbo'
EXTRA="WARNING: slimmc-turbo is memory greedy version of slimmc. Allocated memory is never freed."
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\nvar gcc=\"$GCC\"\n" > ./slimmc/version.nim
nim c -d:release --gc:none -o:./bin/slimmc-turbo ./slimmc/slimmc.nim
rm ./slimmc/version.nim


