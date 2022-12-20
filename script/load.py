import shutil
import argparse
from script.settings import *

def write_fasta_from_json(input_json):
	# Creates a fasta file from json file.
	if os.path.isfile(os.path.join(data_path, "proteindb.fsa")):
		logger.info("Database already exists.")
		return
	else:
		if not os.path.exists(db_path):
			logger.info("Making database folder at {}".format(PATH))
			os.system('mkdir ' + db_path)
		if not os.path.exists(path):
			logger.info("Making database folder {}".format(path))
			os.system('mkdir ' + path)
		if not os.path.exists(data_path):
			logger.info("Making datafolder {}".format(data_path))
			os.system('mkdir ' + data_path)

		with open(os.path.join(PATH, input_json), 'r') as jfile:
			data = json.load(jfile)
			with open(os.path.join(data_path, "proteindb.fsa"), 'a') as fout:
				for i in range(0, len(data)):
					header = (data[i]['dbAMP_ID'])
					sequence = data[i]['Seq']
					fout.write(">{}\n".format(header))
					fout.write("{}\n".format(sequence))
		logger.info("Database extract finish: {}".format(os.path.join(data_path, "proteindb.fsa")))


def main(args):
	if args.debug:
		logger.setLevel(10)

	write_fasta_from_json(args.input_json)


def create_parser():
	parser = argparse.ArgumentParser(prog="macpfinder load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('-i', '--input_json', required=True, help="input file must be in json format! e.g ampfinder_210729.json")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()