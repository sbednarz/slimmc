GTAG=`git describe --abbrev=4 --dirty --always --tags`
GHASH=`git describe --always --dirty`
BUILD=`date`
printf "var gtag=\"$GTAG\"\nvar ghash=\"$GHASH\"\nvar build=\"$BUILD\"\n" > version.nim
nim c -d:release slimmc.nim
rm version.nim

