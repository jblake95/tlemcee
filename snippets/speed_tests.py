from skyfield.api import load, utc
from skyfield.sgp4lib import EarthSatellite
from datetime import datetime
from astropy.table import Table
import numpy as np
import time

from st.site import Site

EPHEM_PATH = '/Users/jblake95/GitHub/tlemcee/ephem_astra_1m_20180618_LaPalma.csv'
TS = load.timescale()

ephem = Table.read(EPHEM_PATH)

line1 = '1 33436U 08057A   18170.03521762  .00000104  00000-0  00000+0 0  9996'
line2 = '2 33436   0.0136 355.3865 0004990 118.0808 185.6030  1.00271245 35197'

obj = EarthSatellite(line1, line2)
site = Site('LaPalma')

# test 1 - without arrays
times = []
for t in ephem['UTC']:
    times.append(datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f'))
times = np.array(times)
print('Computing for {} positions'.format(len(times)))

from astropy.coordinates import SkyCoord
from astropy import units as u

t0 = time.time()
ra, dec = np.zeros((2, len(times)))
for t, dt in enumerate(times):
    r, d, _ = ((obj - site.topocentric).at(TS.utc(dt.replace(tzinfo=utc)))).radec()
    coord = SkyCoord(ra = r.hours,
                    dec = d.degrees,
                    unit=(u.hourangle, u.deg),
                    frame='icrs')
    ra[t] = coord.ra.deg
    dec[t] = coord.dec.deg
t1 = time.time()
print('Test 1 took {} seconds'.format(t1 - t0))

# test 2 - with arrays
times_replaced = []
for t in times:
    times_replaced.append(t.replace(tzinfo=utc))

t0 = time.time()
r, d, _ = ((obj - site.topocentric).at(TS.utc(times_replaced))).radec()
r = r._degrees
d = d.degrees
t1 = time.time()
print('Test 2 took {} seconds'.format(t1 - t0))

# test 3 - precompute time attributes
times_replaced = TS.utc(times_replaced)
times_replaced.MT
times_replaced.gast

t0 = time.time()
r, d, _ = ((obj - site.topocentric).at(times_replaced)).radec()
r = r._degrees
d = d.degrees
t1 = time.time()
print('Test 3 took {} seconds'.format(t1 - t0))

from st.tle import TLE
tle = TLE(line1, line2)
tle.parse_propagation_info(times, 'LaPalma')

t0 = time.time()
ra, dec = tle.propagate_radec()
t1 = time.time()
print('Test 4 took {} seconds'.format(t1 - t0))
