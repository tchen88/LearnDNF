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
for correlation in correlations:
	for diff_method in diff_methods:
		for hypothesis in initial_hypotheses:
			for error in errors:
				for nt in num_terms:
					for nv in zip(num_vars, num_vars_per_term):
						for index in range(5):
							print './LearnDNF.out -s {0:d} -a {1} -m {2} -h {3} -e {4} -t {5} -v {6} -z {7} > outputs/{0:d}.txt'.format(seed_num,
								correlation, diff_method, hypothesis, error, nt, nv[0], nv[1])
							seed_num += 1