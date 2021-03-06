GTAG=`git describe --abbrev=4 --dirty --always --tags`
GHASH=`git rev-parse HEAD`
BUILD=`date`
SYS=`uname`
NIM=`nim --version | head -n 1`
GCC=`gcc --version | head -n 1`

echo '>>> Building slimmc'
PRG='slimmc'
EXTRA=''
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\nvar gcc=\"$GCC\"\n" > version.nim
nim c -d:release -o:slimmc slimmc.nim
rm version.nim

echo ">>> Building slimmc-turbo"
PRG='slimmc-turbo'
EXTRA="*** Warning! Memory greedy version ***"
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\nvar gcc=\"$GCC\"\n" > version.nim
nim c -d:release --gc:none -o:slimmc-turbo slimmc.nim
rm version.nim
