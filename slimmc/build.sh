
GTAG=`git describe --abbrev=4 --dirty --always --tags`
GHASH=`git describe --always --dirty`
BUILD=`date`
SYS=`uname`
NIM=`nim --version | head -n 1`

printf "var gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\nvar sys=\"$SYS\"\nvar nimv=\"$NIM\"\n" > version.nim
nim c -d:release slimmc.nim
rm version.nim

