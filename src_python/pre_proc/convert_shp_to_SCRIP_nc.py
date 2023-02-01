#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : convert_shp_to_SCRIP_nc
Project Name: MPR
Description : converts shapefile format (.shp) to SCRIP netCDF format (.nc) accepted by MPR and vice versa
Author      : Robert Schweppe (robert.schweppe@ufz.de)
Created     : 2019-09-05 11:46
"""

import argparse
import math
import warnings
from itertools import product

# IMPORTS
import geopandas as gpd
import numpy as np
import pandas as pd
import tqdm
import xarray as xr
from shapely.geometry import Polygon

# GLOBAL VARIABLES
MISSING_VALUE = -9999.0
DEFAULT_FULL_INFORMATION = False
DEFAULT_UNITS = 'degrees'
LOOKUP_UNITS = {
    'degrees': 'degrees',
    'radians': 'radians',
}


# FUNCTIONS

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''
        Converting a ESRI polygon shapefile (*.shp) to a SCRIP-formatted netCDF file (*.nc) and vice versa.

        author: Robert Schweppe
        created: Sep 2019''')
    parser.add_argument('-i', '--input_file', action='store',
                        dest='input_file', metavar='input_file',
                        help="path to input file")
    parser.add_argument('-o', '--output_file', action='store',
                        dest='output_file',
                        help="path to output file")
    parser.add_argument('-n', '--name', action='store',
                        dest='name', required=False,
                        help="name for Dimension - "
                             "set this when writing to nc")
    parser.add_argument('-f', '--full_information', action='store_true',
                        default=DEFAULT_FULL_INFORMATION, dest='full_information',
                        help="store all attributes in netcdf file? (Default: {}) -"
                             "set this when writing to nc".format(DEFAULT_FULL_INFORMATION))
    parser.add_argument('-u', '--units', action='store',
                        default=DEFAULT_UNITS, dest='units',
                        help="units of the coordinate system? (Default: {}) - "
                             "set this when writing to nc".format(DEFAULT_UNITS))

    return parser.parse_args()


# CLASSES
class MyGeoDataFrame(gpd.GeoDataFrame):
    def get_SCRIP_vars(self):
        # self.check_type()
        print('getting the centroid lon values')
        # centroid_lon = self.check_longitude(self.get_centroid('x'))
        centroid_lon = self.get_centroid('x')
        print('getting the centroid lat values')
        centroid_lat = self.get_centroid('y')
        corner_lon, corner_lat = self.get_corners_lon()
        return centroid_lon, centroid_lat, corner_lon, corner_lat

    def check_type(self):
        if not self.geom_type == 'Polygon':
            raise

    def get_centroid(self, attr):
        centroid = np.array([getattr(item, attr) for item in self.geometry.centroid])
        return centroid

    def get_corners_lon(self):
        """
        get all the corners of each polygon as a 2d array, if polygon is of type MultiPolygon, only use first polygon
        """
        # first get the number of corners for each polygon
        # subtract 1 because last vertex is double
        print('getting number of vertices for each cell')
        lengths = [len(item.exterior.coords.xy[0]) - 1 if item.geom_type == 'Polygon' else len(
            item[0].exterior.coords.xy[0]) - 1 for item in self.geometry]
        max_length = max(lengths)
        # init the final arrays and set the default missing value
        corner_lon = np.zeros((len(self.geometry), max_length))
        corner_lon[corner_lon == 0] = MISSING_VALUE
        corner_lat = np.zeros((len(self.geometry), max_length))
        corner_lat[corner_lat == 0] = MISSING_VALUE
        # now loop over each polygon and iteratively set the corner values in the target array
        # TODO: sort the nodes so the centroid is always left of the node path
        # https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        for i_item, item in enumerate(self.geometry):
            if item.geom_type == 'Polygon':
                corner_lon[i_item, :lengths[i_item]], corner_lat[i_item, :lengths[i_item]] = \
                    self._check_order(*item.exterior.coords.xy)
            elif item.geom_type == 'MultiPolygon':
                corner_lon[i_item, :lengths[i_item]], corner_lat[i_item, :lengths[i_item]] = \
                    self._check_order(*item[0].exterior.coords.xy)
            else:
                raise
        # corner_lon = self.check_longitude(corner_lon)
        return corner_lon, corner_lat

    def check_longitude(self, lon_arg):
        if lon_arg.min() <= 0:
            warnings.warn('Longitude values are ranging from -180 to 180, converting to range 0 to 360')
            lon_arg = np.where(lon_arg < 0, lon_arg + 360, lon_arg)
        return lon_arg

    def _check_order(self, lons, lats):
        """
        check the order of the polygon vertices and ensure they are counter-clockwise

        Parameters
        ----------
        lons: x-coordinates
        lats: y-coordinates

        Returns
        -------
        lons: checked and counter-clockwise x-coordinates
        lats: checked and counter-clockwise y-coordinates

        """
        # drop the last vertex (is same as first)
        lons = lons[:-1]
        lats = lats[:-1]

        min_index = np.argmin(lats)
        while True:
            # all neighboring indices
            indices = [(min_index - 1) % len(lats), min_index, (min_index + 1) % len(lats)]
            # calculate determinant as here: https://en.wikipedia.org/wiki/Curve_orientation
            det = \
                (lons[indices[1]] - lons[indices[0]]) * (lats[indices[2]] - lats[indices[0]]) - \
                (lons[indices[2]] - lons[indices[0]]) * (lats[indices[1]] - lats[indices[0]])
            if det == 0:
                # the three points are colinear, check next vertices
                min_index = indices[2]
                continue
            elif det > 0:
                # already counter-clockwise
                return lons, lats
            else:
                # sort ascending
                return lons[::-1], lats[::-1]


def handle_latlon_format(ds):
    coord_var = ['south_north', 'west_east']
    # select all data_vars that have dimension grid_size
    shp_vars = [data_var for data_var in ds.data_vars if all([dim in ds[data_var].dims for dim in coord_var])]
    # empty container
    dfs = []
    for data_var in shp_vars:
        # get all the dimension names for this data_var except grid_size
        remaining_dims = [dim for dim in ds[data_var].dims if dim not in coord_var]
        # reshape array, so x-axis is grid_size and remaining dimensions are pandas multiIndex
        if remaining_dims:
            df = ds[data_var].stack(new_index=coord_var).to_pandas().T
            # index are flat index again
            df.index = ['_'.join((str(item) for item in _)) for _ in df.index]
            df.columns = ['{}_{}'.format(data_var, i) for i in df.columns]
        else:
            df = pd.DataFrame(ds[data_var].to_pandas())
            df.columns = [data_var]
        dfs.append(df)
    # we turn the pandas DataFrame in a Geopandas GeoDataFrame and still need geometry
    gdf = gpd.GeoDataFrame(pd.concat(dfs, axis=1), geometry=[[]] * len(dfs[0]))
    # convert units to degree
    # loop over polygons and calculate number of edges and set geometry item
    for i, j in product(ds[coord_var[0]], ds[coord_var[1]]):
        xvals = ds['XLONG_bnds'].isel(**{coord_var[1]: j}).values
        yvals = ds['XLAT_bnds'].isel(**{coord_var[0]: i}).values
        gdf.loc['{}_{}'.format(float(i), float(j)), 'geometry'] = Polygon([
            (xvals[0], yvals[0]),
            (xvals[1], yvals[0]),
            (xvals[1], yvals[1]),
            (xvals[0], yvals[1])
        ])

    return gdf


# SCRIPT
if __name__ == '__main__':
    """
    sys.argv = ['',
                #'-i', '/Users/ottor/nc/Home/projects/2018_MPR_paper/data/04_remap_SCRIP/08_Florida/shapes/WBDHU6_FLORIDA.shp',
                '-i', '/Users/ottor/nc/Home/local_libs/fortran/MPR/src/tests/files/multiple_polygons.shp',
                #'-i', '/Users/ottor/eve/ottor/lib/MPR/src/tests/files/polygons.shp',
                #'-i', '/Users/ottor/nc/Home/local_libs/fortran/MPR/default_output.nc',
                #'-i', '/Users/ottor/nc/Home/local_libs/fortran/MPR/deps/SCRIP/SCRIP/grids/remap_grid_T42.nc',
                #'-i', '/Users/ottor/Downloads/icon_grid_0019_R02B05_G.nc',
                #'-i', '/Users/ottor/nc/Home/projects/2018_MPR_paper/data/04_remap_SCRIP/06_latlon/parameters/'
                #      'cosby_CONUS_soilgrids_default_deg.nc',
                '-o', '/Users/ottor/nc/Home/local_libs/fortran/MPR/src/tests/files/multiple_polygons.nc',
                #'-o', '/Users/ottor/nc/Home/projects/2018_MPR_paper/data/04_remap_SCRIP/02_HUC6/WBDHU6_CONUS_east.nc',
                #'-o', '/Users/ottor/nc/Home/projects/2018_MPR_paper/data/04_remap_SCRIP/08_Florida/shapes/WBDHU6_FLORIDA.nc',
                #'-o', '/Users/ottor/eve/ottor/lib/MPR/src/tests/files/polygons.nc',
                #'-o', '/Users/ottor/nc/Home/projects/2018_MPR_paper/data/04_remap_SCRIP/04_ICON/icon_grid_0019_R02B05_G_CONUS_east.nc',
                #'-n', 'WBDHU6_CONUS',
                '-n', 'polygons',
                '-u', 'degrees']
    """
    # parse the arguments
    args = parse_args()

    if args.input_file.endswith('.shp') and args.output_file.endswith('.nc'):
        # read the shapefile
        mygdf = MyGeoDataFrame(gpd.read_file(args.input_file))
        try:
            units = LOOKUP_UNITS[mygdf.crs['units']]
        except (TypeError, ValueError):
            units = args.units
        centroid_lon, centroid_lat, corner_lon, corner_lat = mygdf.get_SCRIP_vars()

        ds = xr.Dataset(
            data_vars={'grid_corner_lon': (['grid_size', 'grid_corners'], corner_lon),
                       'grid_corner_lat': (['grid_size', 'grid_corners'], corner_lat),
                       'grid_center_lon': (['grid_size'], centroid_lon),
                       'grid_center_lat': (['grid_size'], centroid_lat),
                       },
            attrs={'title': args.name, 'units': units}
        )
        # add units globally and to each variable
        for data_var in ds.data_vars:
            ds[data_var].attrs['units'] = units
            if data_var.startswith('grid_center_'):
                ds[data_var].attrs['bounds'] = data_var.replace('center', 'corner')
            if data_var.startswith('grid_corner_'):
                ds[data_var].attrs['_FillValue'] = MISSING_VALUE

        # add dimension grid_rank and add variable grid_dims(grid_rank)
        # they need to be set to 1 for unstructured grids, as there are not dimensions for polygons
        ds['grid_dims'] = xr.DataArray(np.array([len(centroid_lat)], dtype=int), dims=['grid_rank'])

        # add variable grid_imask(grid_size)
        ds['grid_imask'] = xr.DataArray(np.ones_like(centroid_lon, dtype=int),
                                        dims=['grid_size'],
                                        attrs={'units': 'unitless'},
                                        )
        print('Writing', len(centroid_lat), 'polygons to file', args.output_file, 'with name', args.name)
        ds.to_netcdf(args.output_file)

    elif args.input_file.endswith('.nc') and args.output_file.endswith('.shp'):
        # all SCRIP variable names
        EXCLUDE_DIMS = ['grid_center_lon', 'grid_center_lat', 'grid_corner_lon', 'grid_corner_lat', 'grid_imask']
        # all MPI-ICON based variable names and its SCRIP counterparts
        RENAME_VAR_DICT = {'clon': 'grid_center_lon', 'clat': 'grid_center_lat',
                           'clon_vertices': 'grid_corner_lon', 'clat_vertices': 'grid_corner_lat'}
        # all MPI-ICON based dimension names and its SCRIP counterparts
        RENAME_DIM_DICT = {'cell': 'grid_size', 'nv': 'grid_corners'}
        # all MPI-ICON based data variable names
        SELECT_VARS = ['cell_area_p', 'cell_elevation', 'cell_sea_land_mask']
        # read the netcdf file
        ds = xr.open_dataset(args.input_file)
        # SCRIP format
        mpi_icon_format = all([dim in ds for dim in RENAME_VAR_DICT.keys()])
        scrip_format = all([dim in ds for dim in EXCLUDE_DIMS[:4]])
        if not (mpi_icon_format or scrip_format):
            gdf = handle_latlon_format(ds)
        else:
            # do we have special MPI-ICON names?
            if mpi_icon_format:
                # select and rename
                ds = ds[list(RENAME_VAR_DICT.keys()) + SELECT_VARS].rename({**RENAME_VAR_DICT, **RENAME_DIM_DICT})
                # set the grid_imask property based on sea_land_mask or to default 1
            if 'cell_sea_land_mask' in ds:
                ds['grid_imask'] = xr.where(ds['cell_sea_land_mask'] > 1, 1, 0)
            else:
                ds['grid_imask'] = (('grid_size',), np.ones_like(ds['grid_center_lon'].values, dtype=int))
            # select all data_vars that have dimension grid_size
            shp_vars = [data_var for data_var in ds.data_vars if
                        'grid_size' in ds[data_var].dims and data_var not in EXCLUDE_DIMS]
            if not shp_vars:
                ds['id'] = (('grid_size',), np.array(range(1, len(ds['grid_size']) + 1)))
                shp_vars = ['id']
            # empty container
            dfs = []
            for data_var in shp_vars:
                # get all the dimension names for this data_var except grid_size
                remaining_dims = [dim for dim in ds[data_var].dims if dim != 'grid_size']
                # reshape array, so x-axis is grid_size and remaining dimensions are pandas multiIndex
                if remaining_dims:
                    df = ds[data_var].stack(new_index=remaining_dims).to_pandas()
                    # columns are flat index again
                    df.columns = ['_'.join((str(item) for item in _)) for _ in df.columns]
                else:
                    df = pd.DataFrame(ds[data_var].to_pandas())
                    df.columns = [data_var]
                dfs.append(df)
            # we turn the pandas DataFrame in a Geopandas GeoDataFrame and still need geometry
            gdf = gpd.GeoDataFrame(pd.concat(dfs, axis=1, keys=shp_vars))
            # merge the MultiIndex into a flat index
            gdf.columns = ['_'.join((str(item) for item in _)) for _ in gdf.columns]

            # convert units to degree
            n_corners = len(ds['grid_corners'])
            for var_name in EXCLUDE_DIMS:
                units = ds[var_name].attrs.get('units')
                if units and units.startswith('radian'):
                    ds[var_name] *= 180 / math.pi
            # check for common number of edges for cells
            if (ds['grid_corner_lon'] == -9999).any():
                # loop over polygons and calculate number of edges and set geometry item
                for i_polygon in tqdm.trange(len(ds['grid_size'])):
                    n = n_corners - (ds['grid_corner_lon'].isel(grid_size=i_polygon).values == -9999).sum()
                    gdf.loc[i_polygon, 'geometry'] = Polygon(zip(
                        ds['grid_corner_lon'].isel(grid_size=i_polygon, grid_corners=slice(n)).values,
                        ds['grid_corner_lat'].isel(grid_size=i_polygon, grid_corners=slice(n)).values))
            else:
                # loop over polygons and set geometry item
                for i_polygon in tqdm.range(len(ds['grid_size'])):
                    gdf.loc[i_polygon, 'geometry'] = Polygon(zip(
                        ds['grid_corner_lon'].isel(grid_size=i_polygon).values,
                        ds['grid_corner_lat'].isel(grid_size=i_polygon).values))

        # set WGS84 by default
        gdf.crs = {'init': 'epsg:4326'}
        gdf.to_file(args.output_file)
    else:
        raise Exception('Did not get proper filenames. Script only works for conversion of *.shp->*.nc or vice versa')
