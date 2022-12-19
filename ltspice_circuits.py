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
        elif self.analysis == 'ac_get_data':
            return l.get_data(node), l.get_time()

def run_ac_analysis():
    #netlist = ['Evaluation circuit','V1 1 0 SINE(0 10 1000)', 'R1 1 2 1000', 'C1 2 0 470p', '.tran 10m 1', '.end']
    netlist = [ 'Evaluation circuit', 'Q1 Out 2 3 0 NPN', 'R1 1 Out 6800', 'R2 3 0 1000', 'V1 1 0 9', 
                'V2 5 0 SINE(0 0.05 300)', 'R3 1 2 4700', 'R4 2 0 1000', 
                'C1 3 0 20u', 'C2 2 5 1u', '.model NPN NPN', '.model PNP PNP',
                '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt',
                '.tran 1m 500m',
                '.backanno',
                '.end']

    my_circuit = Circuit(analysis='ac_get_data')
    my_circuit.generate_cir(netlist)
    my_circuit.run_cir()
    out, t = my_circuit.read_output('V(Out)')
    print(t[-250:])
    plt.plot(out[-250:])
    t = np.arange(0, 0.5)
    print(t[:10])
    ref = 4.5*np.sin(2*np.pi*300*t[-250:] + np.pi) + 4.5
    plt.plot(ref)
    plt.show()

def run_dc_mapping():
    netlist = [ 'DC Mapping', 'R1 N001 N002 3300', 'R2 N001 Out 3300', 'R3 N006 N008 150', 
                'R4 N005 N008 150', 'R5 N008 N007 3300', 'Q1 N002 N003 N006 0 2N3904',
                'Q2 Out N004 N005 0 2N3904', 'V1 N003 0 0', 'V2 N004 0 1', 'V3 N001 0 1',
                'V4 N007 0 -1', '.model NPN NPN', '.model PNP PNP', 
                '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt',
                '.dc V1 -3 3 0.01', '.backanno', '.end']
    my_circuit = Circuit(analysis='get_data')
    my_circuit.generate_cir(netlist)
    my_circuit.run_cir()
    out = my_circuit.read_output('V(Out)')

    v_in = np.arange(-3, 3, 0.01)
    v_out = 1/(1 + np.exp(-6*v_in))

    plt.plot(out)
    plt.plot(v_out)

    plt.show()


if __name__ == '__main__':
    import matplotlib.pyplot as plt 
    import numpy as np 
    
    #run_ac_analysis()
    run_dc_mapping()