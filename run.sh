#!/usr/bin/bash
echo 'Warning: This could take up to two days to run.'
echo 'This requires Python and g++ to run. It also requires ~125MB of disk space'
./WriteShellScript.py > LearnDNFScript.sh
./LearnDNFScript.sh
./ScrapeResults.py > Results.csv