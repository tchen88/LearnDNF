#!/usr/bin/env python

correlations = ['None', '0.1', '0.3', '0.5', '0.7', '0.9']
diff_methods = ['Difference', 'Ratio']
initial_hypotheses = ['1', 'Probability ratio', 'Number of terms']
errors = ['0.02', '0.05', '0.1']
num_terms = ['4', '8']
num_vars = ['8', '10']
num_vars_per_term = ['3', '4']
seed_num = 0
print 'Alpha\tDifference method\tInitial hypothesis\tError tolerance\tNumber of terms\tNumber of variables\tVariables per term\tRunning time'
for correlation in correlations:
	for diff_method in diff_methods:
		for hypothesis in initial_hypotheses:
			for error in errors:
				for nt in num_terms:
					for nv in zip(num_vars, num_vars_per_term):
						for index in range(5):
							infile = open('outputs/{0:d}.txt'.format(seed_num))
							for line in infile:
								if line[0:12] == 'Running time':
									line = line.strip('s\n').split(' ')
									print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.3f}'.format(correlation, diff_method, hypothesis, error, nt, nv[0], nv[1], float(line[2]))
							seed_num += 1