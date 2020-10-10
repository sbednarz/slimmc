GTAG=`git describe --abbrev=4 --dirty --always --tags`
GHASH=`git rev-parse HEAD`
BUILD=`date`
SYS=`uname`
NIM=`nim --version | head -n 1`

echo '>>> Building slimmc'
PRG='slimmc'
EXTRA=''
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\n" > version.nim
nim c -d:release -o:slimmc11 slimmc.nim
rm version.nim

echo ">>> Building slimmc-turbo"
PRG='slimmc-turbo'
EXTRA='*** Warning! Memory greedy version ***'
printf "var prg=\"$PRG\"\nvar extra=\"$EXTRA\"\nvar gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\n" > version.nim
nim c -d:release --gc:none -o:slimmc-turbo11 slimmc.nim
rm version.nim
