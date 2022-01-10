#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : split_domain
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 11.08.20 09:23
"""

# IMPORTS
from pyresample.geometry import create_area_def
from pyresample.utils import load_cf_area
from pyresample import AreaDefinition
from pyresample import geo_filter
import pathlib as pl
from typing import Optional, Tuple, Union, Iterable, Generator
import xarray as xr
import numpy as np
from pandas import to_numeric
import argparse
from textwrap import dedent
import multiprocessing as mp
from math import ceil

# GLOBAL VARIABLES
MASK_VARIABLE_NAME = 'ids'
DEFAULT_COMPRESSION = {'zlib': True, 'complevel': 5}
DEFAULT_TARGET_DIR = '~/temp'
DEFAULT_FILE_PATTERN = 'subdomain_{i_grid}.nc'
N_CORE = mp.cpu_count()


# FUNCTIONS
def read_netcdf_file(filename: Union[str, pl.Path], variable: Optional[str] = None) \
        -> Union[xr.Dataset, xr.DataArray]:
    with xr.open_dataset(filename) as data:
        if variable is not None and variable in data.data_vars:
            data = data[variable]
    return data


def read_area_definition(filename: Union[str, pl.Path], **kwargs) \
        -> AreaDefinition:
    """
    see https://pyresample.readthedocs.io/en/latest/geometry_utils.html?highlight=create_area_def#loading-from-netcdf-cf
    """
    area, cf_info = load_cf_area(filename, **kwargs)
    return area


def set_area_definition(**kwargs) \
        -> AreaDefinition:
    area = create_area_def(**kwargs)
    # this needs to be set to make the geo_filter stuff work
    if area.lons is None:
        area.lons, area.lats = area.get_lonlats()
    return area


def create_subdomains_from_maskfile(area: AreaDefinition, maskfile: Union[str, pl.Path]) \
        -> Generator[geo_filter.GridFilter, None, None]:
    basin_ids = read_netcdf_file(maskfile, MASK_VARIABLE_NAME).data
    # check if conforms with area
    if basin_ids.shape != area.shape:
        raise ValueError(f'Shapes of mask ({basin_ids.shape}) and data grid ({area.shape}) do not match')
    unique_ids = np.unique(basin_ids[~np.isnan(basin_ids)])
    for unique_id in unique_ids:
        yield geo_filter.GridFilter(area, (basin_ids == unique_id))


def _get_filter_from_block(area: AreaDefinition, sub_width: float, sub_height: float, block: int) \
        -> geo_filter.GridFilter:
    arr = np.zeros(shape=area.shape, dtype=bool)
    low_x = int((block) * sub_width) % area.width
    up_x = int((block + 1) * sub_width - 1) % area.width + 1
    low_y = int(int((block * sub_width) / area.width) * sub_height) % area.height
    up_y = int(int((block * sub_width) / area.width + 1) * sub_height - 1) % area.height + 1
    # print(block, low_x, up_x, low_y, up_y)
    arr[slice(low_y, up_y), slice(low_x, up_x)] = True
    return geo_filter.GridFilter(area, arr)


def create_subdomains_from_dims(area: AreaDefinition, split_dims: Iterable[str], n_subdomains: int) \
        -> Iterable[geo_filter.GridFilter]:
    if split_dims == ('y',):
        n_subdomains = blocks = min(area.height, n_subdomains)
        sub_width = area.width
        sub_height = area.height / n_subdomains
    elif split_dims == ('x',):
        n_subdomains = blocks = min(area.width, n_subdomains)
        sub_width = area.width / n_subdomains
        sub_height = area.height
    elif split_dims == ('x', 'y'):
        n_subdomains = min(area.width * area.height, n_subdomains)
        cells_per_block = np.sqrt((area.width * area.height) / n_subdomains)
        if cells_per_block > area.width:
            sub_width = area.width
            sub_height = area.height / n_subdomains
            blocks = n_subdomains
        elif cells_per_block > area.height:
            sub_width = area.height
            sub_height = area.width / n_subdomains
            blocks = n_subdomains
        elif area.height < area.width:
            n_larger_side = int(area.width / cells_per_block)
            sub_width = area.width / n_larger_side
            n_smaller_side = int(n_subdomains / n_larger_side)
            sub_height = area.height / n_smaller_side
            blocks = n_larger_side * n_smaller_side
        else:
            n_larger_side = int(area.height / cells_per_block)
            sub_height = area.height / n_larger_side
            n_smaller_side = int(n_subdomains / n_larger_side)
            sub_width = area.width / n_smaller_side
            blocks = n_larger_side * n_smaller_side
    else:
        raise Exception('Provided invalid split_dims')

    for block in range(blocks):
        yield _get_filter_from_block(area, sub_width, sub_height, block)


def create_subdomains_from_dxdy(area: AreaDefinition, split_dims: Iterable[str],
                                dxdy: Tuple[Union[int, float], Union[int, float]],
                                ll: Optional[Tuple[Union[int, float], Union[int, float]]] = None) \
        -> Iterable[geo_filter.GridFilter]:

    def _x(index_arg, *args):
        return index_arg, 0

    def _y(index_arg, *args):
        return 0, index_arg

    def _xy(index_arg, subarea_shape):
        return index_arg % subarea_shape[0], index_arg // subarea_shape[0]

    # get the ll corners
    ll = check_ll(ll, dxdy, area.area_extent[:2])
    dxdy = check_dxdy(dxdy, area.resolution)
    if split_dims == ('y',):
        subdomains_per_dim = (1, ceil((area.area_extent[3] - ll[1]) / dxdy[1]))
        dxdy = (area.area_extent[2] - ll[0], dxdy[1])
        index_func = _x
    elif split_dims == ('x',):
        subdomains_per_dim = (ceil((area.area_extent[2] - ll[0]) / dxdy[0]), 1)
        dxdy = (dxdy[0], area.area_extent[3] - ll[1])
        index_func = _y
    elif split_dims == ('x', 'y'):
        subdomains_per_dim = (ceil((area.area_extent[2] - ll[0]) / dxdy[0]),
                              ceil((area.area_extent[3] - ll[1]) / dxdy[1]))
        index_func = _xy
    else:
        raise Exception('Provided invalid split_dims')
    for i in range(subdomains_per_dim[0] * subdomains_per_dim[1]):
        indices = index_func(i, subdomains_per_dim)
        yield _get_filter_from_dxdy(area, ll, dxdy, indices)


def _get_filter_from_dxdy(area: AreaDefinition, ll: Tuple[Union[int, float], Union[int, float]],
                          dxdy: Tuple[Union[int, float], Union[int, float]], i: Tuple[int, int]) \
        -> geo_filter.GridFilter:
    llx, lly, urx, ury = (
        ll[0] + dxdy[0] * i[0],
        ll[1] + dxdy[1] * i[1],
        ll[0] + dxdy[0] * (i[0] + 1),
        ll[1] + dxdy[1] * (i[1] + 1)
    )
    # # this code snippet attempted to use pyresample based implementation, but that did not work
    # proj_def = area.crs.to_dict() if hasattr(area, 'crs') else area.proj_str
    # sub_area = create_area_def('temp', projection=proj_def, area_extent=(llx, lly, urx, ury), resolution=dxdy)
    # sub_area.crs = area.crs
    proj_x, proj_y = area.get_proj_coords()
    x_mask = (proj_x - 0.5 * area.resolution[0] >= llx) & \
             (proj_x + 0.5 * area.resolution[0] <= urx)
    y_mask = (proj_y - 0.5 * area.resolution[1] >= lly) & \
             (proj_y + 0.5 * area.resolution[1] <= ury)
    return geo_filter.GridFilter(area, x_mask & y_mask)


def check_ll(ll: Iterable[Union[int, float]], dxdy: Iterable[Union[int, float]], area_ll: Iterable[Union[int, float]]):
    if ll is None:
        ll = area_ll
    else:
        # if ll outside
        ll = np.array(ll).astype(float)
        dxdy = np.array(dxdy)
        ll += (np.array(area_ll) - ll) // dxdy * dxdy
    return tuple(ll)


def check_dxdy(dxdy: Iterable[Union[int, float]], area_dxdy: Iterable[Union[int, float]]):
    new_dxdy = np.array(dxdy).astype(float)
    for i in range(len(dxdy)):
        # check if dxdy are at least the resolution of the are
        if new_dxdy[i] < area_dxdy[i]:
            new_dxdy[i] = area_dxdy[i]
        # other than that, make them a multiple of the resolution (nearest lesser integer)
        if new_dxdy[i] % area_dxdy[i]:
            new_dxdy[i] = new_dxdy[i] - new_dxdy[i] % area_dxdy[i]
    return tuple(new_dxdy)


def split_subdomains(areas: Iterable[geo_filter.GridFilter],
                     area: AreaDefinition,
                     data: Union[xr.Dataset, xr.DataArray],
                     dims: Optional[Tuple[str, str]] = None,
                     file_pattern=DEFAULT_FILE_PATTERN,
                     nprocs=1,
                     use_mask=True) \
        -> Generator[Tuple[Union[xr.Dataset, xr.DataArray], str], None, None]:
    # use arg, fall back to data.dims
    dims = dims or data.dims
    # create coordinates from dims, if not yet existing
    proj_vectors = area.get_proj_vectors()[::-1]
    for idim, dim in enumerate(dims):
        if dim not in data.coords and len(data[dim]) == len(proj_vectors[idim]):
            data.coords[dim] = proj_vectors[idim]
    # loop over subdomains, simply indexing them with a counter
    for i_filter, grid_filter in enumerate(areas, 1):
        # # the pyresample produces a 1D array of only the valid data, lons and lats
        # area_filtered, data_filtered = grid_filter.filter(self.area, self.data.data)
        # xr.where with drop=True requires a boolean DataArray instead of np.array
        subdomain = data.where(xr.DataArray(grid_filter.get_valid_index(area),
                                            dims=dims),
                               drop=True)
        mask = subdomain.isnull()
        if use_mask and mask.all():
            subdomain = None
        elif use_mask and mask.any():
            for dim in dims:
                subdomain = subdomain.dropna(dim, how='all')
        filename = file_pattern.format(i_grid=i_filter)
        yield subdomain, filename


def write_netcdf_file(obj: Union[xr.Dataset, xr.DataArray], filename: Union[str, pl.Path], compress=True):
    encoding = {}
    if compress:
        compress_kwargs = {}
        # set compression
        if isinstance(obj, xr.Dataset):
            compress_kwargs = {var: DEFAULT_COMPRESSION for var in obj.data_vars.keys()}
        elif hasattr(obj, 'name'):
            compress_kwargs = {obj.name: DEFAULT_COMPRESSION}
        # make more sophisticated if required
        encoding = compress_kwargs
    obj.to_netcdf(filename, encoding=encoding)


def parse_args(args: Optional[str] = None):
    """parse the arguments passed from the console and display help with '-h' switch"""
    parser = argparse.ArgumentParser(
        description=dedent('''\
        Create spatial subdomains for netcdf files
        
        author: Robert Schweppe
        created: Aug 2020'''))
    parser.add_argument('--filename', dest='filename', type=str,
                        help="path of file to be split into subdomains")
    parser.add_argument('--grid_kwargs', dest='grid_kwargs',
                        action=StoreDictKeyPair,
                        nargs="+", default=None, help=dedent('''\
    if grid cannot be inferred from file, 
    create it using https://pyresample.readthedocs.io/en/latest/geometry_utils.html?highlight=create_area_def#loading-from-netcdf-cf
    provide the according key-value pairs like key1=value1 key2=value2 ...'''))
    parser.add_argument('--target_dir', dest='target_dir', default=None, help="folder to write output files to")
    parser.add_argument('--maskfile', dest='maskfile', default=None, help=dedent('''\
    path of file containing masks for subdomains
    - needs to be in valid netCDF4 format and needs to contain variable 'ids' with the subdomain identifiers
    - flag missing values accordingly with '_FillValue'
    - the shape of the variable needs to be in accordance with the Area Definition of the data variable'''))
    parser.add_argument('--y', dest='y', default=None, help="name of y-dimension")
    parser.add_argument('--x', dest='x', default=None, help="name of x-dimension")
    parser.add_argument('--n_subdomains', dest='n_subdomains', default=None, type=int, help="number of subdomains")
    parser.add_argument('--dxdy', dest='dxdy', default=None, type=lambda x: to_numeric(x.split(',')),
                        help="list of subdomain width and subdomain height, comma-separated")
    parser.add_argument('--ll', dest='ll', default=None, type=lambda x: to_numeric(x.split(',')),
                        help="list of subdomain x left boundary, y lower boundary, comma-separated")
    parser.add_argument('--split_dims', dest='split_dims', default=None, type=lambda arg: tuple(arg.split(',')),
                        help=dedent("""\
    dimensions to be split into subdomains, pass one of 'x', 'y', 'x,y'
    - x -> split along x dimension (into columns)  
    - y -> split along y dimension (into rows)
    - x,y -> split into blocks, you might get less n_subdomains than specified  
    """))
    parser.add_argument('--variable', dest='variable', default=None,
                        help="variable of netcdf file (file_name) to be selected")
    parser.add_argument('--trim_to_mask', dest='trim_to_mask', action='store_true',
                        help="trim subdomain tiles individually to their closest bounding boxes (dropna is applied)")

    return parser.parse_args(args)


# CLASSES
class StoreDictKeyPair(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super(StoreDictKeyPair, self).__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        my_dict = {}
        # loop over nargs
        for kv in values:
            k, v = kv.split("=", 1)
            if ',' in v:
                v = to_numeric(v.split(','))
            my_dict[k] = v
        setattr(namespace, self.dest, my_dict)


class SubdomainCreater():

    def __init__(self, filename, variable=None, target_dir=DEFAULT_TARGET_DIR, grid_kwargs=None,
                 split_dims=None, n_subdomains=None, maskfile=None, subdomain_areas=None, y=None, x=None,
                 dxdy=None, ll=None, nprocs=N_CORE):
        # load the data contents
        filename = pl.Path(filename).expanduser()
        self.data = read_netcdf_file(filename=filename, variable=variable)

        # set the grid information
        self.area = None
        if grid_kwargs is not None:
            self.area = set_area_definition(**grid_kwargs)
        else:
            try:
                self.area = read_area_definition(filename=filename, variable=variable, y=y, x=x)
            # there is a bug in pyresample.py raising an AttributeError
            except (ValueError, AttributeError):
                pass
        if not isinstance(self.area, AreaDefinition):
            raise Exception(f'Unable to set AreaDefinition, pass valid grid_kwargs')

        self.dims = None
        if y is not None and x is not None:
            self.dims = (y, x)

        # init the subdomains
        self._init_subdomains(split_dims, n_subdomains, maskfile, subdomain_areas, dxdy, ll)

        # set the target directory
        self.target_dir = target_dir

        self.nprocs = nprocs

    def _init_subdomains(self, split_dims, n_subdomains, maskfile, subdomain_areas, dxdy, ll):
        if split_dims is not None and n_subdomains is not None:
            self.subdomain_areas = create_subdomains_from_dims(self.area, split_dims, n_subdomains)
        elif split_dims is not None and dxdy is not None:
            self.subdomain_areas = create_subdomains_from_dxdy(self.area, split_dims, dxdy, ll)
        elif maskfile is not None:
            self.subdomain_areas = create_subdomains_from_maskfile(self.area, maskfile)
            # split along rows
        elif subdomain_areas is not None:
            self.subdomain_areas = subdomain_areas
        else:
            raise Exception(f'Did not provide valid information on how to split subdomains')

    def split(self, *args, **kwargs):
        for subdomain, filename in split_subdomains(areas=self.subdomain_areas, area=self.area, data=self.data,
                                                    dims=self.dims, nprocs=self.nprocs, *args, **kwargs):
            if subdomain is None:
                print(f'skipping empty subdomain: {filename}')
                continue
            print(f'writing subdomain to filename: {filename}')
            write_netcdf_file(subdomain, pl.Path(self.target_dir, filename).expanduser())


# SCRIPT
if __name__ == '__main__':
    # # example args
    # example_args = '--filename ~/temp/all.nc --variable ids --target_dir ~/temp/out --grid_kwargs area_id=EPSG:4326 projection=EPSG:4326 shape=1800,3600 area_extent=-180,-90,180,90 --maskfile ~/temp/all.nc'
    # example_args = '--filename ~/temp/land_cover.nc --variable land_cover --target_dir ~/temp/out --y lat --x lon --dxdy 40000,40000 --ll 4220000,5320000 --split_dims x,y --grid_kwargs area_id=EPSG:31468 projection=EPSG:31468 shape=1720,1520 area_extent=4226000,5324000,4378000,5496000 --trim_to_mask'
    args = parse_args(
        # example_args.split()
    )
    my_subdomains = SubdomainCreater(
        filename=args.filename,
        variable=args.variable,
        grid_kwargs=args.grid_kwargs,
        target_dir=args.target_dir,
        maskfile=args.maskfile,
        y=args.y,
        x=args.x,
        n_subdomains=args.n_subdomains,
        dxdy=args.dxdy,
        ll=args.ll,
        split_dims=args.split_dims,
    )
    my_subdomains.split(use_mask=args.trim_to_mask)
