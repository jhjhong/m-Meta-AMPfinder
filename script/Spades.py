from script.settings import *
import shutil

class Spades(object):

    def __init__(self, input1, input2=None, output_dir=None, num_threads=16):
        self.input1 = input1
        self.input2 = input2
        self.output_dir = output_dir
        self.num_threads = num_threads

    def __repr__(self):
        return "Metaspades({}".format(self.__dict__)

    def run(self):
        f1_path, f1_name = os.path.split(self.input1)
        q1_name = "{}.temp.highqual.fq".format(f1_name.split(".")[0])

        stdout = "2> /dev/null"

        if self.input2:
            f2_path, f2_name = os.path.split(self.input2)
            q2_name = "{}.temp.highqual.fq".format(f2_name.split(".")[0])

            cmd = "metaspades.py -1 {q1_file} -2 {q2_file} -o {out} -t {num_threads} {stdout}" \
            .format(
                q1_file=os.path.join(f1_path, q1_name),
                q2_file=os.path.join(f2_path, q2_name),
                out=os.path.join(self.output_dir, "temp.assembly"),
                num_threads=self.num_threads,
                stdout=stdout
            )
        else:
            cmd = "megahit --presets meta-large -r {q1_file} -o {out} --min-contig-len 1000 -t {num_threads} {stdout}" \
                .format(
                    q1_file=os.path.join(f1_path, q1_name),
                    out=os.path.join(self.output_dir, "temp.assembly"),
                    num_threads=self.num_threads,
                    stdout=stdout
                )
        # print(cmd)
        os.system(cmd)
        
        # move final file
        original = os.path.join(self.output_dir, "temp.assembly/final.contigs.fa")
        target = os.path.join(self.output_dir, "final.contigs.fasta")
        shutil.copyfile(original, target)