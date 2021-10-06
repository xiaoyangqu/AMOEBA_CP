#! bin/bash/

# usage: sh calc_amoeba_energy.sh {txyz-file} {key-file}

export OPENMMHOME=/data/home/qxy/tinker/bin

function ana()
{
	txyz=$1
	key=$2
	$OPENMMHOME/analyze ${txyz} -k ${key} d > ${txyz%.*}.log 2>err.log &
}

ana $1 $2
