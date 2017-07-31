#!/usr/bin/env python

correlations = ['None', '0.1', '0.3', '0.5', '0.7', '0.9']
diff_methods = ['Difference', 'Ratio']
initial_hypotheses = ['1', 'p/q', '2t - 1']
errors = ['0.02', '0.05', '0.1']
num_terms = ['4', '8']
num_vars = ['8', '10']
num_vars_per_term = ['3', '4']
print ('DNF expression\tp\tNumber of terms\tNumber of variables\tVariables per term\tCorrelation\tDifference method\t'
	'Initial hypothesis\tError tolerance\tRunning time (s)\tNumber of iterations\tNumber of coefficients\tTimed out?')
for nt in num_terms:
	for nv in zip(num_vars, num_vars_per_term):
		for correlation in correlations:
			for error in errors:
				for diff_method in diff_methods:
					for hypothesis, hypothesis_abbrev in zip(initial_hypotheses, ['1', 'p', 't']):
						for index in range(5):
							infile = open('outputs/{0}T{1}V{2}C{3}M{4}H{5:d}_{6}e.txt'.format(
								nt, nv[0], 'n' if correlation == 'None' else correlation, diff_method[0].lower(), hypothesis_abbrev, index, error))
							initial_prob = None
							num_iter = 0
							num_coef = 0
							is_timeout = False
							expression = None
							for line in infile:
								if initial_prob == None and line[0:13] == 'Current error':
									line = line.strip('\n').split(' ')
									initial_prob = 1 - float(line[2])
								elif line[0:15] == 'Loop iteration:':
									num_iter += 1
								elif line[0:16] == 'New coefficient!':
									num_coef += 1
								elif line[0:10] == 'Timed out:':
									is_timeout = True
								elif line[0:4] == 'DNF:':
									line = line.strip('s\n').split(':')
									expression = line[1].strip()
								elif line[0:13] == 'Running time:':
									line = line.strip('s\n').split(' ')
									print '{0}\t{1:.4f}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9:.3f}\t{10}\t{11}\t{12}'.format(
										expression, initial_prob, nt, nv[0], nv[1], correlation, diff_method, hypothesis, error, float(line[2]),
										num_iter, num_coef, 'Yes' if is_timeout else 'No')