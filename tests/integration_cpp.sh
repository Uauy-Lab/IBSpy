#/bin/sh

lines="jagger arinalrfor julius lancer landmark mace norin61 spelta stanley sy_mattis"

data_path="./tests/data"
kmer_db=$data_path/test4B.jagger.kmerGWAS_k31
expected_output="tests/data/ibspy_integration_test.txt"
output="tests/ibspy_out.txt"

for l in $lines; do
	ref="$data_path/test4B.$l.fa"
	echo $ref
	./IBScpp/build/IBScpp -d $kmer_db -r $ref -k 31 -w 3000 -p 2
done > $output

def_ret=`diff $output $expected_output`

echo "Differences: $def_ret"
echo "Code: ${#def_ret}"
exit ${#def_ret}

