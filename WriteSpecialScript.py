#!/usr/bin/env python

correlation = 'n'
diff_methods = ['d', 'r']
initial_hypotheses = ['1', 't']
errors = ['0.001', '0.02']
nt = '8'
nv = ['8', '3']
seed_num = 1000
print '#!/usr/bin/bash'
print 'g++ LearnDNF.cpp -o LearnDNF.out'
for diff_method in diff_methods:
	for hypothesis in initial_hypotheses:
		for index in range(100):
			for error in errors:
				print './LearnDNF.out -s {0:d} -a {1} -m {2} -h {3} -e {4} -t {5} -v {6} -z {7} > outputs/special{5}T{6}V{1}C{2}M{3}H{8}_{4}e.txt'.format(seed_num,
					correlation, diff_method, hypothesis, error, nt, nv[0], nv[1], index)
			seed_num += 1