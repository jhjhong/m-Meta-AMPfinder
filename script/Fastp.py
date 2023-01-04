from script.settings import *


class Fastp(object):

    def __init__(self, short1, short2=None, output_dir=None, num_threads=16):
        self.short1 = short1
        self.short2 = short2

        self.output_dir = output_dir
        self.num_threads = num_threads

    def __repr__(self):
        return "Fastp({}".format(self.__dict__)

    def run(self):
        f1_path, f1_name = os.path.split(self.short1)
        q1_name = "{}.temp.highqual.fsa".format(f1_name.split(".")[0])
        json= "{}.temp.fastp.json".format(f1_name.split(".")[0])
        html= "{}.temp.fastp.html".format(f1_name.split(".")[0])

        if self.short2:
            f2_path, f2_name = os.path.split(self.short2)
            q2_name = "{}.temp.highqual.fsa".format(f2_name.split(".")[0])
            short2 = "-I {short2} -O {out2}".format(short2=self.short2, out2=os.path.join(self.output_dir, q2_name))
        
        stdout = "2> /dev/null"

        cmd = "fastp -i {short1} -o {out1} {short2} -n 0 -w {num_threads} -j {json} -h {html} {stdout}" \
            .format(
                short1=self.short1,
                out1=os.path.join(self.output_dir, q1_name),
                short2=short2,
                num_threads=self.num_threads,
                json=os.path.join(self.output_dir, json),
                html=os.path.join(self.output_dir, html),
                stdout=stdout
            )
        # print(cmd)
        os.system(cmd)