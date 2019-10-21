# calcualte the normalized z-score for each of the top 10 models using single-linkage clustering

import argparse
import numpy as np
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import linkage, fcluster

# parse the arguments from command line
parser = argparse.ArgumentParser(description = 'Calculate the normalized z-score for top 10 models from VESPER. Normalized z-scores for top 10 models are written into the output file.')
parser.add_argument('-i', '--input', required = True, action = 'store', dest = 'input_file', help = 'Required. Name of input file.')
parser.add_argument('-c', type = float, action = 'store', dest = 'cutoff', default = 0.2, help = 'Optional. Clustering cutoff ranging from 0 to 1. Default = 0.2.')
parser.add_argument('-o', '--output', action = 'store', dest = 'out_name', help = 'Optional. Name of output file. If not specified, the output file would be named as input filename followed by .normzscore.')

args = parser.parse_args()

input_file = filename = args.input_file
cutoff = args.cutoff

if not args.out_name:
	output_name = input_file + '.normzscore'
else:
	output_name = args.out_name

# read in score information in input file
with open(filename) as zscore_file:
	model_line_start = ['#0', '#1', '#2', '#3', '#4', '#5', '#6', '#7', '#8', '#9']
	model_score_info = np.empty(10)

	zscore_content = zscore_file.read()
	if 'Score=' in zscore_content:
		score_list = []
		top_model_score_list = []
		for line in zscore_content.split('\n'):
			if line[0:2] == 'R ':
				score = float(line.split()[-1])
				score_list.append(score)
			elif line[0:2] in model_line_start:
				model_num = int(line[1])
				score = float(line.split()[-3])
				model_score_info[model_num] = score

		if np.any(np.isnan(score_list)) or np.any(np.isnan(model_score_info)):
			print(score_list)
			print(model_score_info)
			print('Error: NaN value exists in ' + input_file)
			print('Calculation of normalized z-score is not performed')
			print('Check ' + input_file + ' to make sure VESPER has run properly')

		else:
			score_list = np.array(score_list)
			max_d = cutoff * (max(score_list) - min(score_list))
			single_score_list = [[i] for i in score_list]

			# single linkage clustering
			clusters = fcluster(linkage(single_score_list, 'single'), max_d, criterion = 'distance')

			# calculate mean and std of the largest cluster
			cluster_count = Counter(clusters)
			first_cluster_index = sorted(cluster_count, key = cluster_count.get, reverse = True)[0]
			first_cluster = score_list[clusters == first_cluster_index]
			adjust_mean, adjust_std = np.mean(first_cluster), np.std(first_cluster)

			output = open(output_name, 'w')
			output.write('Normalized z-score for top 10 models:\n')
			print('Normalized z-score for top 10 models:')
			# calculate normalized z-score for each model
			for i in np.arange(10):
				old_score = model_score_info[i]
				norm_zscore = (old_score - adjust_mean)/adjust_std

				output.write('#' + str(i) + '\t' + str(norm_zscore) + '\n')
				print('#' + str(i) + '\t' + str(norm_zscore))
		
			output.close()
		
	else:
		print('Error: There is no DOT score information in ' + input_file)
		print('Check ' + input_file + ' to make sure VESPER has run properly')


