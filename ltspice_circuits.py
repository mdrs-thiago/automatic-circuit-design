import ltspice 
import os 
import subprocess


class Circuit:

    def __init__(self, analysis, filename:str = 'circuit', spice_path:str = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe'):
        self.analysis =  analysis
        self.filename = filename
        self.circuit_path = os.path.join('circuits',f'{filename}.cir')
        self.result_path = os.path.join('circuits',f'{filename}.raw')
        self.spice_path = spice_path

    def generate_cir(self, components: list):
        netlist = ' \n '.join(components)

        with open(self.circuit_path, 'w') as f:
            f.write(netlist)

    def run_cir(self):
        subprocess.call(f"{self.spice_path} --Run -b {self.circuit_path}");

    def read_output(self, node):
        l = ltspice.Ltspice(self.result_path)
        l.parse()
        if self.analysis == 'get_data':
            return l.get_data(node)
