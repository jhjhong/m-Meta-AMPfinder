from script.settings import *


class Megahit(object):

    def __init__(self, input1, input2=None, output_dir=None, num_threads=16):
        self.input1 = input1
        self.input2 = input2
        self.output_dir = output_dir
        self.num_threads = num_threads

    def __repr__(self):
        return "Megahit({}".format(self.__dict__)

    def run(self):
        f1_path, f1_name = os.path.split(self.input1)
        q1_name = "{}.temp.highqual.fsa".format(f1_name.split(".")[0])

        stdout = "2> /dev/null"

        if self.input2:
            f2_path, f2_name = os.path.split(self.input2)
            q2_name = "{}.temp.highqual.fsa".format(f2_name.split(".")[0])

            cmd = "megahit --presets meta-large -1 {q1_file} -2 {q2_file} --min-contig-len 1000 -o {out} -t {num_threads} {stdout}" \
            .format(
                q1_file=os.path.join(self.output_dir, q1_name),
                q2_file=os.path.join(self.output_dir, q2_name),
                out=os.path.join(self.output_dir, "temp.assembly"),
                num_threads=self.num_threads,
                stdout=stdout
            )
        else:
            cmd = "megahit --presets meta-large -r {q1_file} -o {out} --min-contig-len 1000 -t {num_threads} {stdout}" \
                .format(
                    q1_file=os.path.join(self.output_dir, q1_name),
                    out=os.path.join(self.output_dir, "temp.assembly"),
                    num_threads=self.num_threads,
                    stdout=stdout
                )
        print(cmd)
        # os.system(cmd)