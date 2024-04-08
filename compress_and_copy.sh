ls *.dat > ll
l=`cat ll | wc -l`
for i in `seq 1 $l`;
do
  line=`head -n $i ll | tail -n 1`
  zstd -9 "$line" -o "$line".zstd
  sha512sum "$line".zstd > "$line".zstd.sha512sum
done
cp -v *.zstd.sha512sum $1
cp -vf timing* $1
cp -v status* $1
cp -v *.zstd $1
echo "`hostname` has finished writeout"
