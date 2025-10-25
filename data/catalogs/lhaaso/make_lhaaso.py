"""
““” ““” Script to convert dictionary `data` of LHAASO sources to an
ECSV format table.

The twelve LHAASO sources are from:
https://www.nature.com/articles/s41586-021-03498-z ““”

““”

"""



import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from gammapy.modeling.models import PowerLawSpectralModel, LogParabolaSpectralModel
from gammapy.utils.scripts import make_path
from feupy.utils.table import pad_list_to_length
from feupy.utils.units import FRAME_ICRS, UNIT_DEG, DEFAULT_ENERGY_UNIT, DEFAULT_SED_UNIT
from feupy.utils.constants import CU
# from feupy.scripts.ipynb_to_gallery import convert_ipynb_to_gallery

def format_table(table):
    """Format table columns for better readability."""
    for column in table.colnames:
        if column.startswith(("Flux", 'sed_dnde')):
            table[column].format = ".3e"
        elif column.startswith(("e_min", "e_max", "e_ref", "sed_e", "sqrt_ts", "norm", "ts", "stat")):
            table[column].format = ".3f"
    return table


def read_sed_data(file, length):
    """Read spectral energy distribution (SED) data from a file."""
    e_ref_list, dnde_list, dnde_errn, dnde_errp = [], [], [], []
    with open(file, 'r') as f:
        for line in f:
            values = line.strip().split()
            if len(values) >= 2:
                try:
                    e_ref = float(values[0]) / 1e12  # Convert to TeV
                    f_val = float(values[1])
                    err_up = float(values[2]) if len(values) > 2 else f_val
                    err_down = float(values[3]) if len(values) > 3 else err_up
                    if f_val / err_down < 1:
                        continue
                    e_ref_list.append(e_ref)
                    dnde_list.append(f_val)
                    dnde_errn.append(err_down)
                    dnde_errp.append(err_up)
                except ValueError:
                    continue
    return (
        pad_list_to_length(length, e_ref_list),
        pad_list_to_length(length, dnde_list),
        pad_list_to_length(length, dnde_errn),
        pad_list_to_length(length, dnde_errp),
    )

# Constants
BIBCODE = '2021Natur.594...33C'
TAG = 'LHAASO'

REF_URL = "https://doi.org/10.1038/s41586-021-03498-z"
DATA_FILES = {
    'LHAASO J1825-1326': 'J1825_KM2A_201209.dat',
    'LHAASO J1908+0621': 'J1908_KM2A_201209.dat',
    'LHAASO J2226+6057': 'J2228_KM2A_201209.dat',
}
SED_TYPE = 'e2dnde'
DEFAULT_LENGTH = 10
REFERENCE_ENERGY = 100 * u.TeV

CATALOG = {
    'LHAASO J0534+2202': {
        'position': SkyCoord(83.55, 22.05, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 17.8,
        'Energy_max':(0.88, 0.11)*u.Unit('PeV'),
        'Flux_100TeV':(1.00, 0.14)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=1.00*CU,
            reference=100*u.TeV,
        )},
    'LHAASO J1825-1326': {
        'position': SkyCoord(276.45, -13.45, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 16.4,
        'Energy_max':(0.42, 0.16)*u.Unit('PeV'),
        'Flux_100TeV': (3.57, 0.52)*CU,
        'spectral_model': LogParabolaSpectralModel(
            alpha=0.92,
            amplitude='1e-12 cm-2 s-1 TeV-1',
            reference=10*u.TeV,
            beta=1.19)},  
    'LHAASO J1839-0545': {
        'position': SkyCoord(279.95, -5.75, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 7.7,
        'Energy_max': (0.21, 0.05)*u.Unit('PeV'),
        'Flux_100TeV': (0.70, 0.18)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.70*CU,
            reference=100*u.TeV,
        )},  
    'LHAASO J1843-0338': {
        'position': SkyCoord(280.75, -3.65, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 8.5,
        'Energy_max' :(0.26, 0.10)*u.Unit('PeV'),
        'Flux_100TeV': (0.73,0.17)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.73*CU,
            reference=100*u.TeV,
        )}, 
    'LHAASO J1849-0003': {
        'position': SkyCoord(282.35, -0.05, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 10.4,
        'Energy_max': (0.35, 0.07)*u.Unit('PeV'),
        'Flux_100TeV': (0.74, 0.15)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.73*CU,
            reference=100*u.TeV,
        )}, 
    'LHAASO J1908+0621': {
        'position': SkyCoord(287.05, 6.35, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 17.2,
        'Energy_max': (0.44, 0.05)*u.Unit('PeV'),
        'Flux_100TeV': (1.36, 0.18)*CU,
        'spectral_model': LogParabolaSpectralModel(
            alpha=2.27,
            amplitude='1e-12 cm-2 s-1 TeV-1',
            reference=10*u.TeV,
            beta=0.46)},
    'LHAASO J1929+1745': {
        'position': SkyCoord(292.25, 17.75, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 7.4,
        'Energy_max': (0.71, 0.07)*u.Unit('PeV'),
        'Flux_100TeV': (0.38, 0.09)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.38*CU,
            reference=100*u.TeV,
        )},  
    'LHAASO J1956+2845': {
        'position': SkyCoord(299.05, 28.75, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 7.4,
        'Energy_max': (0.42, 0.03)*u.Unit('PeV'),
        'Flux_100TeV': (0.41, 0.09)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.41*CU,
            reference=100*u.TeV,
        )},
    'LHAASO J2018+3651': {
        'position': SkyCoord(304.75, 36.85, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 10.4,
        'Energy_max': (0.27, 0.02)*u.Unit('PeV'),
        'Flux_100TeV': (0.50, 0.10)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.5*CU,
            reference=100*u.TeV,
        )},
    'LHAASO J2032+4102': {
        'position': SkyCoord(308.05, 41.05, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 10.5,
        'Energy_max': (1.42, 0.13)*u.Unit('PeV'),
        'Flux_100TeV': (0.54, 0.10)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.54*CU,
            reference=100*u.TeV,
        )},
    'LHAASO J2108+5157': {
        'position': SkyCoord(317.15, 51.95, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 8.3,
        'Energy_max': (0.43, 0.05)*u.Unit('PeV'),
        'Flux_100TeV': (0.38, 0.09)*CU,
        'spectral_model': PowerLawSpectralModel(            
            amplitude=0.38*CU,
            reference=100*u.TeV,
        )},
    'LHAASO J2226+6057': {
        'position': SkyCoord(336.75, 60.95, unit=UNIT_DEG, frame=FRAME_ICRS),
        'Significance_100TeV': 13.6,
        'Energy_max': (0.57, 0.19)*u.Unit('PeV'),
        'Flux_100TeV': (1.05, 0.16)*CU,
        'spectral_model': LogParabolaSpectralModel(
            alpha=1.56,
            amplitude='1e-12 cm-2 s-1 TeV-1',
            reference=10*u.TeV,
            beta=0.88)}
    }

# Data Preparation


# Generate Table Columns
catalog_table_meta = {'catalog_name': TAG,
                      "SED_TYPE": SED_TYPE,
                      "reference": REF_URL
}
# Prepare data for the table
source_names = []
ra_list = []
dec_list = []
Significance_100TeV_list = []
Energy_max_list = []
Energy_max_error_list = []
Flux_100TeV_list = []
Flux_error_100TeV_list = []
spec_model_type = []


# Iterate over the sources to extract values
for source_name, source_data in CATALOG.items():
    source_names.append(source_name)
    
    ra_list.append(source_data['position'].ra.deg)
    dec_list.append(source_data['position'].dec.deg)
    Significance_100TeV_list.append(source_data['Significance_100TeV'])
        
    Energy_max_list.append(source_data['Energy_max'][0])
    
    Energy_max_error_list.append(source_data['Energy_max'][1].to_value(u.PeV))
    Flux_100TeV_list.append(source_data['Flux_100TeV'][0].to_value(u.Unit('cm-2 s-1 TeV-1')))
    Flux_error_100TeV_list.append(source_data['Flux_100TeV'][1].to_value(u.Unit('cm-2 s-1 TeV-1')))
    
    # Extract spectral model details
    model = source_data['spectral_model']
    if isinstance(model, PowerLawSpectralModel):
        spec_model_type.append('pl')
    elif isinstance(model, LogParabolaSpectralModel):
        spec_model_type.append('lp')
    else:
        spec_model_type.append('Unknown')
        
table = Table([source_names, ra_list, dec_list, Significance_100TeV_list, Energy_max_list, Energy_max_error_list, 
               Flux_100TeV_list, Flux_error_100TeV_list, spec_model_type],
              names=('source_name', 'ra', 'dec', 'Significance_100TeV', 'Energy_max', 'Energy_max_error', 
                     'Flux_100TeV', 'Flux_error_100TeV', 'spec_type'))
table.meta.update(catalog_table_meta)


e_ref_list, dnde_list, dnde_err_list, dnde_errn_list, dnde_errp_list, dnde_ul_list, is_ul_list = [], [], [], [], [], [], []
for _ in table:
    if _['spec_type'] == 'pl':
        e_ref = pad_list_to_length(DEFAULT_LENGTH, [REFERENCE_ENERGY.value])
        dnde = pad_list_to_length(DEFAULT_LENGTH, [_['Flux_100TeV']])*REFERENCE_ENERGY**2
        dnde_err = pad_list_to_length(DEFAULT_LENGTH, [_['Flux_error_100TeV']])*REFERENCE_ENERGY**2
        dnde_errn = pad_list_to_length(DEFAULT_LENGTH, [])
        dnde_errp = pad_list_to_length(DEFAULT_LENGTH, [])
        dnde_ul = pad_list_to_length(DEFAULT_LENGTH, [])
        is_ul = np.full(DEFAULT_LENGTH, None, dtype=bool).tolist()
    else: 
        source_name  = _['source_name']
        file_path = DATA_FILES[source_name]
        e_ref, dnde, dnde_errn, dnde_errp =  read_sed_data(file_path, length=DEFAULT_LENGTH)
        dnde_err = pad_list_to_length(DEFAULT_LENGTH, [])
        dnde_ul = pad_list_to_length(DEFAULT_LENGTH, [])
        is_ul = np.full(DEFAULT_LENGTH, None, dtype=bool).tolist()
        
    e_ref_list.append(e_ref)
    dnde_list.append(dnde)
    dnde_err_list.append(dnde_err)
    dnde_errn_list.append(dnde_errn)
    dnde_errp_list.append(dnde_errp)
    dnde_ul_list.append(dnde_ul)
    is_ul_list.append(is_ul)
    
col1 = Column(
    data = e_ref_list,
    name = 'sed_e_ref',
    description='Reference energy', unit=DEFAULT_ENERGY_UNIT[SED_TYPE])

col2 = Column(
    data = dnde_list,
    name = 'sed_dnde',
    description='Differential flux (dnde) SED value', unit=DEFAULT_SED_UNIT[SED_TYPE])

col3 = Column(
    data = dnde_err_list,
    name = 'sed_dnde_err',
    description='Differential flux (dnde) SED errors', unit=DEFAULT_SED_UNIT[SED_TYPE])

col4 = Column(
    data = dnde_errn_list,
    name = 'sed_dnde_errn',
    description='Differential flux (dnde) SED negative errors', unit=DEFAULT_SED_UNIT[SED_TYPE])

col5 = Column(
    data = dnde_errp_list,
    name = 'sed_dnde_errp',
    description='Differential flux (dnde) SED positive errors', unit=DEFAULT_SED_UNIT[SED_TYPE])

col6 = Column(
    data = dnde_ul_list,
    name = 'sed_dnde_ul',
    description='Differential flux (dnde) SED upper limit', unit=DEFAULT_SED_UNIT[SED_TYPE])

col7 = Column(
    data = is_ul_list,
    name = 'sed_is_ul',
    description='Whether data is an upper limit.', unit=  u.Unit(""))
table.add_columns([col1, col2, col3, col4, col5, col6, col7])
table['spec_reference'] = 100*u.TeV

table = table[
    'source_name',
    'ra',
    'dec',
    'Significance_100TeV',
    'Energy_max',
    'Energy_max_error',
    'Flux_100TeV',
    'Flux_error_100TeV',
    'sed_e_ref',
    'sed_dnde',
    'sed_dnde_err',
    'sed_dnde_errn',
    'sed_dnde_errp',
    'sed_dnde_ul',
    'sed_is_ul',
    'spec_type',
    'spec_reference'
]


# Format Table
table = format_table(table)

new_names =['source_name',
 'ra',
 'dec',
 'Significance_100TeV',
 'Energy_max',
 'Energy_max_error',
 'Flux_100TeV',
 'Flux_error_100TeV',
 'sed_e_ref',
 'sed_e2dnde',
 'sed_e2dnde_err',
 'sed_e2dnde_errn',
 'sed_e2dnde_errp',
 'sed_e2dnde_ul',
 'sed_is_ul',
 'spec_type',
 'spec_reference']
names = table.colnames

table.rename_columns(names, new_names)

# Write Table
table.write('lhaaso_catalog.ecsv', format='ascii.ecsv', overwrite=True)

# def format_info_basic_catalog_table(table):
#     keys = table.colnames
    
#     table["source_name"].description = "Source Name"
#     table["ra"].description = "Right Ascension in the icrs frame"
#     table["ra"].format = ".3f"
#     table["dec"].description = "Declination in the icrs frame"
#     table["dec"].format = ".3f"
#     for key in keys:
#         if key.startswith(("Flux")) or key == 'amplitude':
#             table[key].format = ".3e"
            
#     return table





