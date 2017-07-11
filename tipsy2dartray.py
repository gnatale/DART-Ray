import pynbody
import h5py

sim = input('Please input the name of your simulation file: ')

s = pynbody.load(sim)
s.physical_units()

with h5py.File("output.hdf5", 'w') as f:
    f.create_dataset('starcoord', data=s.s['pos'])
    f.create_dataset('gascoord', data=s.g['pos'])
    f.create_dataset('mstar', data=s.s['mass'])
    f.create_dataset('mgas', data=s.g['mass'])
    f.create_dataset('agestar', data=s.s['age'])
    f.create_dataset('gastemp', data=s.g['temp'])
    f.create_dataset('fehstar', data=s.s['feh'])
    f.create_dataset('fehgas', data=s.g['feh'])
    f.create_dataset('ofegas', data=s.g['ofe'])
