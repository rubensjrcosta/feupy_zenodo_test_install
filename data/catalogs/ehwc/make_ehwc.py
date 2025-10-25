"""
This script converts dictionry `data` to ECSV format table.

The HAWC sources in this source catalog are from
https://doi.org/10.1103/PhysRevLett.124.021102

"""
import numpy as np
import astropy.units as u

from astropy.table import Table, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from gammapy.utils.scripts import make_path

from feupy.utils.table import pad_list_to_length


reference =  'https://doi.org/10.1103/PhysRevLett.124.021102'
sed_type = "e2dnde"

data = {
    "eHWC J0534+220": {
        "ra": 83.61 * u.deg, "ra_err": 0.02 * u.deg,
        "dec": 22.00 * u.deg, "dec_err": 0.03 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": np.nan,  # Point Source
        "flux_56TeV": 1.2 * 1e-14 * (u.ph / (u.cm**2 * u.s)), 
        "flux_56TeV_err": 0.2 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 12.0,
        "sqrt_TS_100TeV": 4.44,
        "nearest_2HWC_source": "J0534+220", 
        "distance_to_2HWC": 0.02 * u.deg,
    },
    "eHWC J1909-193": {
        "ra": 272.46 * u.deg, 
        "ra_err": 0.13 * u.deg,
        "dec": -19.34 * u.deg, 
        "dec_err": 0.14 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.34 * u.deg,
        "flux_56TeV": 2.4 * 1e-14 * (u.ph / (u.cm**2 * u.s)), 
        "flux_56TeV_err": 0.5 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 6.97,
        "sqrt_TS_100TeV": 4.82,
        "nearest_2HWC_source": "J1809-190", 
        "distance_to_2HWC": 0.03 * u.deg,
    },
    "eHWC J1825-134": {
        "ra": 276.40 * u.deg, 
        "ra_err": 0.06 * u.deg,
        "dec": -13.37 * u.deg, 
        "dec_err": 0.10 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.35 * u.deg,
        "flux_56TeV": 4.6 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "flux_56TeV_err": 0.5 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 7.33,
        "sqrt_TS_100TeV": 7.30,
        "nearest_2HWC_source": "J1825-134", 
        "distance_to_2HWC": 0.07 * u.deg,
    },
    "eHWC J1839-057": {
        "ra": 279.77 * u.deg, "ra_err": 0.12 * u.deg,
        "dec": -5.71 * u.deg, "dec_err": 0.11 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.34 * u.deg,
        "flux_56TeV": 1.5 * 1e-14 * (u.ph / (u.cm**2 * u.s)), "flux_56TeV_err": 0.3 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 7.03,
        "sqrt_TS_100TeV": 6.06,
        "nearest_2HWC_source": "J1837-065", "distance_to_2HWC": 0.96 * u.deg,
    },
    "eHWC J1842-035": {
        "ra": 280.72 * u.deg, "ra_err": 0.15 * u.deg,
        "dec": -3.15 * u.deg, "dec_err": 0.15 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.31 * u.deg,
        "flux_56TeV": 1.1 * 1e-14 * (u.ph / (u.cm**2 * u.s)), "flux_56TeV_err": 0.3 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 6.63,
        "sqrt_TS_100TeV": 2.70,
        "nearest_2HWC_source": "J1849-001", "distance_to_2HWC": 0.44 * u.deg,
    },
    "eHWC J1850+001": {
        "ra": 282.59 * u.deg, "ra_err": 0.21 * u.deg,
        "dec": 0.14 * u.deg, "dec_err": 0.12 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.37 * u.deg,
        "flux_56TeV": 1.3 * 1e-14 * (u.ph / (u.cm**2 * u.s)), 
        "flux_56TeV_err": 0.2 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 6.63,
        "sqrt_TS_100TeV": 3.04,
        "nearest_2HWC_source": "J1849-001", "distance_to_2HWC": 0.20 * u.deg,
    },
    "eHWC J1907+063": {
        "ra": 286.91 * u.deg, "ra_err": 0.10 * u.deg,
        "dec": 6.30 * u.deg, "dec_err": 0.12 * u.deg,
        "spectral_model_type": "lp",
        "extension_56TeV": 0.52 * u.deg,
        "flux_56TeV": 2.8 * 1e-14 * (u.ph / (u.cm**2 * u.s)), 
        "flux_56TeV_err": 0.4 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 8.04,
        "sqrt_TS_100TeV": 7.30,
        "nearest_2HWC_source": "J1908+063", "distance_to_2HWC": 0.16 * u.deg,
    },
    "eHWC J2019+368": {
        "ra": 304.95 * u.deg, "ra_err": 0.07 * u.deg,
        "dec": 36.36 * u.deg, "dec_err": 0.08 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.38 * u.deg,
        "flux_56TeV": 3.6 * 1e-14 * (u.ph / (u.cm**2 * u.s)), "flux_56TeV_err": 0.4 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 9.63,
        "sqrt_TS_100TeV": 4.85,
        "nearest_2HWC_source": "J2019+367", "distance_to_2HWC": 0.02 * u.deg,
    },
    "eHWC J2030+412": {
        "ra": 307.74 * u.deg, "ra_err": 0.09 * u.deg,
        "dec": 41.23 * u.deg, "dec_err": 0.07 * u.deg,
        "spectral_model_type": "ecpl",
        "extension_56TeV": 0.18 * u.deg,
        "flux_56TeV": 0.9 * 1e-14 * (u.ph / (u.cm**2 * u.s)), "flux_56TeV_err": 0.2 * 1e-14 * (u.ph / (u.cm**2 * u.s)),
        "sqrt_TS_56TeV": 6.43,
        "sqrt_TS_100TeV": 3.07,
        "nearest_2HWC_source": "J2031+415", "distance_to_2HWC": 0.34 * u.deg,
    },

}


# Step 2: Prepare basic columns for Astropy Table
source_names = list(data.keys())
ra = [data[source]["ra"] for source in source_names]
ra_err = [data[source]["ra_err"] for source in source_names]
dec = [data[source]["dec"] for source in source_names]
dec_err = [data[source]["dec_err"] for source in source_names]
spec_type = [data[source]["spectral_model_type"] for source in source_names]

flux_56TeV = [data[source]["flux_56TeV"] for source in source_names]
flux_56TeV_err = [data[source]["flux_56TeV_err"] for source in source_names]
sqrt_TS_56TeV = [data[source]["sqrt_TS_56TeV"] for source in source_names]
sqrt_TS_100TeV = [data[source]["sqrt_TS_100TeV"] for source in source_names]
nearest_2HWC_source = [data[source]["nearest_2HWC_source"] for source in source_names]
distance_to_2HWC = [data[source]["distance_to_2HWC"] for source in source_names]

table = Table(
    [
        source_names,
        ra,
        ra_err,
        dec,
        dec_err, 
        spec_type,
        flux_56TeV,
        flux_56TeV_err,
        sqrt_TS_56TeV,
        sqrt_TS_100TeV,
        nearest_2HWC_source,
        distance_to_2HWC
    ],
    names=[
        "source_name",
        "ra",
        "ra_err",
        "dec",
        "dec_err",
            "spec_type",
        "Flux_56TeV",
        "Flux_56TeV_err",
        "sqrt_TS_56TeV",
        "sqrt_TS_100 TeV",
        "nearest_2HWC_source",
        "distance_to_2HWC"
    ]
)
table['Flux_56TeV'].format = ".3e"
table.meta['catalog_name'] = 'eHWC'
table.meta["SED_TYPE"] = sed_type
table.meta['comments'] = [f'reference: {reference}']

table

# Step 3: Add flux points for specific sources

# Arrays for flux point values
sed_e_ref, sed_e2dnde, sed_e2dnde_errp, sed_e2dnde_errn, sed_e2dnde_ul, sed_is_ul = [], [], [], [], [], []

np_nan_list = pad_list_to_length(10, [])


for source_name in source_names:
    if source_name in ["eHWC J1825-134", "eHWC J1907+063", "eHWC J2019+368"]:
        file_path = f'$PYTHONPATH/data/catalogs/ehwc/{source_name.replace(" ", "_")}.fits'
        table_fits = Table.read(make_path(file_path))
        sed_e_ref.append(table_fits["e_ref"])
        sed_e2dnde.append(table_fits["e2dnde"])
        sed_e2dnde_errp.append(table_fits["e2dnde_errp"])
        sed_e2dnde_errn.append(table_fits["e2dnde_errn"])
        sed_e2dnde_ul.append(table_fits["e2dnde_ul"])
        sed_is_ul.append(table_fits["is_ul"])
    else:
        sed_e_ref.append((np_nan_list))
        sed_e2dnde.append(np_nan_list)
        sed_e2dnde_errp.append(np_nan_list)
        sed_e2dnde_errn.append(np_nan_list)
        sed_e2dnde_ul.append(np_nan_list)
        sed_is_ul.append(np_nan_list)

# Step 4: Add flux point columns to the table
table["sed_e_ref"] = Column(data=sed_e_ref, unit="TeV", description="Reference Energy", format=".3f")
table["sed_e2dnde"] = Column(data=sed_e2dnde, unit="TeV cm^-2 s^-1", description="SED value", format=".3e")
table["sed_e2dnde_errp"] = Column(data=sed_e2dnde_errp, unit="TeV cm^-2 s^-1", description="SED positive error", format=".3e")
table["sed_e2dnde_errn"] = Column(data=sed_e2dnde_errn, unit="TeV cm^-2 s^-1", description="SED negative error", format=".3e")
table["sed_e2dnde_ul"] = Column(data=sed_e2dnde_ul, unit="TeV cm^-2 s^-1", description="SED upper limit", format=".3e")
table["sed_is_ul"] = Column(data=sed_is_ul,description="SED upper limit flag")

# # Step 5: Save or print the final table
table.write(f'ehwc_catalog.ecsv', format='ascii.ecsv', overwrite=True)
table





models.write(make_path("$PYTHONPATH/data/hawc/2020PhRvL.124b1102A/models.yaml"))

# from feupy.scripts.ipynb_to_gallery import convert_ipynb_to_gallery

# convert_ipynb_to_gallery('Untitled.ipynb', "make_psrcat.py")


# !pyflakes make_psrcat.py