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
        startupinfo = subprocess.STARTUPINFO()
        # startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        # startupinfo.wShowWindow = subprocess.SW_HIDE
        self.startupinfo = startupinfo

    def generate_cir(self, components: list):
        netlist = ' \n '.join(components)

        with open(self.circuit_path, 'w') as f:
            f.write(netlist)

    def run_cir(self):
        subprocess.call(f"{self.spice_path} --Run -b {self.circuit_path}", startupinfo=self.startupinfo)

    def read_output(self, node):
        l = ltspice.Ltspice(self.result_path)
        l.parse()
        if self.analysis == 'get_data':
            return l.get_data(node)

if __name__ == '__main__':
    import matplotlib.pyplot as plt 
    netlist = ['Evaluation circuit','V1 1 0 SINE(0 10 1000)', 'R1 1 2 1000', 'C1 2 0 470p', '.tran 10m 1', '.end']

    my_circuit = Circuit(analysis='get_data')
    my_circuit.generate_cir(netlist)
    my_circuit.run_cir()
    out = my_circuit.read_output('V(2)')
    plt.plot(out[-250:])
    plt.show()