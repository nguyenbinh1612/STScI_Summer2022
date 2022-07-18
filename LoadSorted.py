import foggie.utils.foggie_load as fog
import yt
from yt import YTQuantity,YTArray
import numpy as np
import h5py

halo = 'Tempest'
hnum = '008508'

#Mogget
snap_name = '/Users/bnguyen/Downloads/DD2427/DD2427'
halo_c_v_name = '/Users/bnguyen/FOGGIE/foggie/halo_infos/'+hnum+'/nref11c_nref9f/halo_c_v'
trackname = '/Users/bnguyen/FOGGIE/foggie/halo_tracks/'+hnum+'/nref11n_selfshield_15/halo_track_200kpc_nref9'
masses_dir = '/Users/bnguyen/FOGGIE/foggie/halo_infos/'+hnum+'/nref11c_nref9f/'
hsfile = '/Users/bnguyen/Downloads/'+halo+'_RD0042_allhalostardata.h5'

f = h5py.File(hsfile,'r')
pids = f['particle_IDs'][:]
hids = f['host_IDs'].asstr()[:]
f.close()

def StarParts(pfilter, data):
    return data[("all", "particle_type")] == 2
    
add_particle_filter("stars", function=StarParts, filtered_type='all', requires=["particle_type"])

ds, region = fog.foggie_load(snap_name, trackname, find_halo_center=True, halo_c_v_name=halo_c_v_name, disk_relative=True, \
                            particle_type_for_angmom='young_stars', do_filter_particles=True, gravity=False,\
                            region='refine_box',masses_dir=masses_dir)
                            
ds.add_particle_filter('stars')
ads = ds.all_data()

allstarmasses = ads['stars','particle_mass']

# Create version of star allstarmasses that contains only halo stars and
# lists the masses in the same order as the halostardata file
allstarinds = ads['stars','particle_index']
index = np.argsort(allstarinds) # put particle IDs for all stars in numerical order
sorted_allstars = allstarinds[index]
sorted_index = np.searchsorted(sorted_allstars,pids)
# Create a mask that will reorder an array so that it's in the same order as the halostardata
# file and will include only halo stars
pindex = np.take(index,sorted_index,mode="clip")
mask = allstarinds[pindex] != pids
res = np.ma.array(pindex,mask=mask)
massarr = allstarmasses[np.ma.compressed(res)].in_units('Msun')
# masarr now lists the masses of the stars in the same order that the halostardata file lists their
# particle IDs, host IDs, creation times, etc
