#!/usr/bin/bash
echo 'This requires Python and g++ to run.'
./WriteSpecialScript.py > SpecialDNFScript.sh
./SpecialDNFScript.sh
./ScrapeSpecialResults.py > SpecialResults.csv