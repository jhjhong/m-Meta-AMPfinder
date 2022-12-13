import shutil
import argparse
from script.settings import *

def write_fasta_from_json(args):
	# working_directory = os.getcwd()
	# if os.path.isfile(os.path.join(self.db, "proteindb.fsa")):
	# 		# logger.info("Database already exists.")
	# 		return
	# else:
		working_directory = os.getcwd()
		with open(os.path.join(working_directory, args.input_json), 'r') as jfile:
			data = json.load(jfile)
			with open(os.path.join(working_directory, "database.fasta"), 'a') as fout:
				for i in range(0, len(data)):
					header = (data[i]['dbAMP_ID'])
					sequence = data[i]['Seq']
					fout.write(">{}\n".format(header))
					fout.write("{}\n".format(sequence))
		print("make db finish")

def main(args):
	write_fasta_from_json(args)


def create_parser():
	parser = argparse.ArgumentParser(prog="macpfinder load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('-i', '--input_json', required=True, help="input file must be in json format! e.g ampfinder_210729.json")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()