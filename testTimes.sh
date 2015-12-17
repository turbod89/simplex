#! /bin/sh

init=8;
fin=10;

res="times.csv";
echo "" > $res;


#	echo "Creating matrices";
#	for i in $(seq $init $fin)
#	do
#	  a=$(echo "2^$i" | bc);
#	  ./spmatrix $a $a -m -1 -M 1 -c 4 | ./spconvert > $(echo "examples/matrix$a.mat");
#	done

echo "Computing time";
for i in $(seq $init $fin)
do
  a=$(echo "2^$i" | bc);
  f=$(echo "examples/matrix$a.mat");

  # t=$( ( cat $f $f | time ./spmultiply > /dev/null ) 2>&1 > /dev/null | egrep "user" | sed -e "s/\([0-9\.]*\)user\(.*\)/\1/g" );
  # t=$( ( cat $f | time ./spldu > /dev/null ) 2>&1 > /dev/null | egrep "user" | sed -e "s/\([0-9\.]*\)user\(.*\)/\1/g" );
  t=$( ( cat $f | time ./spldu > /dev/null ) 2>&1 > /dev/null | egrep "user" | sed -e "s/\([0-9\.]*\)user\(.*\)/\1/g" );

  echo "$a;$t;" >> $res;
  echo "$a $t";

done

exit;

for i in $(seq $init $fin)
do
  a=$(echo "2^$i" | bc);
  rm $(echo "examples/matrix$a.mat");
done
