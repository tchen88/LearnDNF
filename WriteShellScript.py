#!/usr/bin/env python

correlations = ['n', '0.1', '0.3', '0.5', '0.7', '0.9']
diff_methods = ['d', 'r']
initial_hypotheses = ['1', 'p', 't']
errors = ['0.02', '0.05', '0.1']
num_terms = ['4', '8']
num_vars = ['8', '10']
num_vars_per_term = ['3', '4']
seed_num = 0
print '#!/usr/bin/bash'
print 'g++ LearnDNF.cpp -o LearnDNF.out'
for nt in num_terms:
	for nv in zip(num_vars, num_vars_per_term):
		for correlation in correlations:
			for diff_method in diff_methods:
				for hypothesis in initial_hypotheses:
					for index in range(5):
						for error in errors:
							print './LearnDNF.out -s {0:d} -a {1} -m {2} -h {3} -e {4} -t {5} -v {6} -z {7} > outputs/{5}T{6}V{1}C{2}M{3}H{8}_{4}e.txt'.format(seed_num,
								correlation, diff_method, hypothesis, error, nt, nv[0], nv[1], index)
						seed_num += 1