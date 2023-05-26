import argparse
import glob
import shutil
from script.settings import *

# clean other files left over
def clean():
	files = glob.glob(os.path.join(path,"*"))
	for f in files:
		remove_directory(f)
		if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["adraft","xml","fsa","draft","pyc","log"]:
			os.remove(f)

		if os.path.isdir(f) == False:
			if os.path.isfile(f) == True and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["py","md"]:
				pass
			else:
				if os.path.isfile(f):
					logger.info("Remove: {}".format(f))
					os.remove(f)

    # clean data files
	data_files = glob.glob(os.path.join(data_path,"*"))
	for datafile in data_files:
		if os.path.isfile(datafile) and os.path.basename(datafile) not in [".gitignore","__init__.py"]:
			logger.info("Remove: {}".format(datafile))
			os.remove(datafile)
	logger.info("Cleaned directory: {}".format(data_path))

    # clean db files
	if os.path.exists(db_path):
		logger.info("Remove: {}".format(db_path))
		shutil.rmtree(db_path)
	logger.info("Cleaned directory: {}".format(db_path))

def remove_directory(directory_path):
	if os.path.isdir(directory_path):
		try:
			logger.info("Remove dir: {}".format(directory_path))
			shutil.rmtree(directory_path)
		except OSError as e:
			logger.error("{} - {}" % (e.filename, e.strerror))

#remove temporary file
def main(args):
	if args.debug:
		logger.setLevel(10)
	
	clean()


def create_parser():
	parser = argparse.ArgumentParser(prog="mampfinder clean", description="{} - {} - Clean".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
