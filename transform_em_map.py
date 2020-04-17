# transform target map according to the rotation and translation of each of the top 10 models in VESPER output

import os
import argparse

def write_chimera_command(command_file, map1, map2, r_vector, t_vector, output_name):
	command_file = open(command_file, 'w')
	command_file.write('from chimera import runCommand as rc\n\n')
	command_file.write('rc("open ' + map1 + '")\n')
	command_file.write('rc("open ' + map2 + '")\n')

	# rotate map
	command_file.write('rc("turn x ' + r_vector[0] + ' center 0,0,0 coord #0 models #1")\n')
	command_file.write('rc("turn y ' + r_vector[1] + ' center 0,0,0 coord #0 models #1")\n')
	command_file.write('rc("turn z ' + r_vector[2] + ' center 0,0,0 coord #0 models #1")\n')

	# translate map
	command_file.write('rc("move ' + ','.join(t_vector) + ' coord #0 model #1")\n')

	# save transformed map
	command_file.write('rc("vop #1 resample onGrid #0")\n')
	command_file.write('rc("volume #2 save ' + output_name + '")\n')
	command_file.close()


# parse the arguments from command line
parser = argparse.ArgumentParser(description = 'Transform the target map given rotation and translation information in the VESPER output file.')
parser.add_argument('-i1', '--input1', required = True, action = 'store', dest = 'ref_map', help = 'Required. Name of the reference map file.')
parser.add_argument('-i2', '--input2', required = True, action = 'store', dest = 'target_map', help = 'Required. Name of the target map file.')
parser.add_argument('-t', required = True, action = 'store', dest = 'vesper_result', help = 'Required. Name of the result file from VESPER.')
parser.add_argument('-odir', action = 'store', dest = 'out_dir', help = 'Optional. Directory for the transformed target map files. If not specified, the transformed target map files would be written to the current directory')

args = parser.parse_args()

ref_map, target_map = args.ref_map, args.target_map
vesper_result = args.vesper_result

if not args.out_dir:
	out_dir = './'
else:
	out_dir = args.out_dir

command_file = 'chimera_command.py'

# extract rotation and translation information from vesper_result
with open(vesper_result) as result:
	model_line_start = ['#0', '#1', '#2', '#3', '#4', '#5', '#6', '#7', '#8', '#9']

	for line in result:
		if line[0:2] in model_line_start:
			model_num = str(int(line[1]) + 1)
			r_info = line.split()[1:4]
			t_info = line.split()[13:16]

			rotation_vector= [r_info[0][3:], r_info[1], r_info[2][:-1]]
			translation_vector = [t_info[0][3:], t_info[1], t_info[2][:-1]]

			# transform target map
			out_name = out_dir + 'target_transform_model_' + model_num + '.mrc'
			write_chimera_command(command_file, ref_map, target_map, rotation_vector, translation_vector, out_name)
			run_command = 'chimera --silent --nogui ' + command_file 
			os.system(run_command)

			# remove chimera command file
			os.system('rm ' + command_file)
			os.system('rm ' + command_file + 'c')


	


