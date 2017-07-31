# Monotone Boolean function learner Version 0.1.1

## Features:
* Possible inputs:
	+ Can take random input with seeds determined by the Seeds.txt file (which can be edited)
	+ Can take monotone DNF expressions as input
* Options for random input:
	+ Number of variables
	+ Number of terms
	+ Number of variables per term
	+ Correlation between terms (optional)
* Options for all input:
	+ Distance function
	+ Initial hypothesis

## Requirements:
* bash
* g++
* Python interpreter
* 125 MB disk space
* Two days of CPU time

## How to use:
* To replicate the results in my paper, simply run ./run.sh
	+ Output will be contained in Results.csv
* To replicate the results on special expressions, run ./runSpecial.sh
* Otherwise, you can manually compile and run LearnDNF.cpp with your own options

## TODO:
* Support truth tables as input
* Add documentation for how the options work
