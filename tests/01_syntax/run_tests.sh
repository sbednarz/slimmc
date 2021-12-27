#!/bin/bash


do_testA () {

	echo -n "01_SYNTAX $desc ..."
	eval "$cmd > $out"

	if [ $? -ne 0 ]
	then
		echo " FAIL"
		echo "Error running command: $cmd > $out. Here is output:"
		cat $out
		exit 1
	else
		c=`$check`
		if [ -z "$c" ]
		then
		      echo " OK"
		else
		      echo " FAIL"
		      echo "Error running command: $cmd > $out."
		      echo "The output is different than expected. Here are details:"
		      echo $c
		      exit 1
		fi
	fi
}

desc='TEST01 Checking if minimal model file is readable'
out='01_minimal.out'
cmd='../../bin/slimmc -t 01_minimal.model'
check='diff 01_minimal.req 01_minimal.out'
do_testA

desc='TEST02 Checking default values'
out='02_default.out'
cmd='../../bin/slimmc -t 02_default.model'
check='diff 02_default.req 02_default.out'
do_testA

desc='TEST03 Checking parameters setting'
out='03_settings.out'
cmd='../../bin/slimmc -t 03_settings.model'
check='diff 03_settings.req 03_settings.out'
do_testA

desc='TEST04 Checking breakpoints setting'
out='04_breakpoints.out'
cmd='../../bin/slimmc -t 04_breakpoints.model'
check='diff 04_breakpoints.req 04_breakpoints.out'
do_testA


