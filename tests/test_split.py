import pytest
import xarray as xr
import numpy as np
import pathlib as pl
from utils.split_domain import SubdomainCreater, DEFAULT_FILE_PATTERN, read_netcdf_file, check_ll
from itertools import product

ARRAY_SHAPE = (9, 18)
DIM_NAMES = ('y', 'x')
RESOLUTIONS = [
    ((15, 180), [(slice(None), slice(j, j + 1)) for j in range(ARRAY_SHAPE[1])]),
    ((40, 180), [(slice(None), slice(j, j + 2)) for j in range(0, ARRAY_SHAPE[1], 2)]),
    ((42, 180), [(slice(None), slice(j, j + 2)) for j in range(0, ARRAY_SHAPE[1], 2)]),
    ((42.0, 180), [(slice(None), slice(j, j + 2)) for j in range(0, ARRAY_SHAPE[1], 2)]),
    ((360, 180), [(slice(None), slice(None))]),
    # now split only y-dim, slices are inverted, as data array stores x ascendingly and y descendingly
    ((360, 25), [(slice(j, j + 1), slice(None)) for j in range(ARRAY_SHAPE[0]-1, -1, -1)]),
    # range for y-dim is slice(5,9), slice(1,5), slice(0,1)
    ((360, 80), [(slice(max(j, 0), j + 4), slice(None)) for j in range(ARRAY_SHAPE[0]-4, -4, -4)]),
    ((40, 80), [(slice(max(i, 0), i + 4), slice(j, j + 2)) for i, j in product(range(ARRAY_SHAPE[0]-4, -4, -4),
                                                                       range(0, ARRAY_SHAPE[1], 2))]),
]
N_SUBDOMAINS_BLOCKS = [
    (2, [(slice(None), slice(0, 9)), (slice(None), slice(9, 18))]),
    (7, [(slice(0, 4), slice(0, 6)), (slice(0, 4), slice(6, 12)), (slice(0, 4), slice(12, 18)),
         (slice(4, 9), slice(0, 6)), (slice(4, 9), slice(6, 12)), (slice(4, 9), slice(12, 18))]),
    (200, [(slice(i, i + 1), slice(j, j + 1)) for i, j in product(range(ARRAY_SHAPE[0]), range(ARRAY_SHAPE[1]))]),
]
N_SUBDOMAINS_COLUMNS = [
    (2, [(slice(None), slice(0, 9)), (slice(None), slice(9, 18))]),
    (7, [(slice(None), slice(0, 2)), (slice(None), slice(2, 5)), (slice(None), slice(5, 7)),
         (slice(None), slice(7, 10)), (slice(None), slice(10, 12)), (slice(None), slice(12, 15)),
         (slice(None), slice(15, 18))]),
    (20, [(slice(None), slice(j, j + 1)) for j in range(ARRAY_SHAPE[1])]),
]
GRID_KWARGS = [
    (4326, [-180, -90, 180, 90]),
    (31468, [4226000, 5324000, 4244000, 5331000]),
]
XARRAY_TYPES = [
    ('DataArray', True),
    ('Dataset', False)
]


@pytest.fixture(
    scope="class",
    params=XARRAY_TYPES)
def xarray_type(request):
    da = xr.DataArray(data=np.fromfunction(lambda y, x: x * y, ARRAY_SHAPE),
                      dims=DIM_NAMES,
                      name='data')
    kwargs = {}
    if request.param[1]:
        kwargs['variable'] = 'data'
    if request.param[0] == 'Dataset':
        yield da.to_dataset(), kwargs
    else:
        yield da, kwargs


@pytest.fixture(
    scope="class",
    params=GRID_KWARGS)
def grid_kwargs(request):
    projection = f'EPSG:{request.param[0]}'
    yield dict(
        area_id=projection,
        projection=projection,
        # rows, cols
        shape=ARRAY_SHAPE,
        # xl, yl, xu, yu
        area_extent=request.param[1],
    )


@pytest.fixture(
    scope="class",
    params=N_SUBDOMAINS_BLOCKS)
def n_subdomains_blocks(request):
    yield request.param


@pytest.fixture(
    scope="class",
    params=N_SUBDOMAINS_COLUMNS)
def n_subdomains_columns(request):
    yield request.param


@pytest.fixture(
    scope="class",
    params=RESOLUTIONS)
def dxdy(request):
    yield request.param


class TestClass:
    def test_integration_mask(self, tmp_path, xarray_type, grid_kwargs):
        """
        - load data from netcdf file (ideally parse coordinate information and return xr.DataSet)
          - use Stephan's use case and Germany dataset
        - put Grid in a specific container (pyresample.AreaDescription)
        - build filter based on
          - a list of bounding boxes
          - a list of masks
          - a number of subdomains and a type how to do so (rows, cols, blocks)
        - split the data into the subdomains
        - put the meta information on there
        """
        data_file = pl.Path(tmp_path, 'data.nc')
        xarray_type[0].to_netcdf(data_file)

        mask = xr.DataArray(data=np.ones(ARRAY_SHAPE),
                            dims=DIM_NAMES,
                            name='ids')
        mask.loc[dict(x=slice(1, 4), y=slice(0, 2))] = 2
        mask.loc[dict(x=slice(4, 6), y=slice(2, 5))] = 2
        mask_file = pl.Path(tmp_path, 'mask.nc')
        mask.to_netcdf(mask_file)

        subs = SubdomainCreater(filename=data_file,
                                grid_kwargs=grid_kwargs,
                                target_dir=tmp_path,
                                maskfile=mask_file,
                                y=DIM_NAMES[0],
                                x=DIM_NAMES[1],
                                **xarray_type[1])

        # make sure missing values are set
        subs.split(use_mask=False)

        for i_grid in range(1, 3):
            sub_path = pl.Path(tmp_path, DEFAULT_FILE_PATTERN.format(i_grid=i_grid))
            assert sub_path.exists()
            sub_path = read_netcdf_file(sub_path, **xarray_type[1])
            assert sub_path.equals(xarray_type[0].where(mask == i_grid, drop=True))

    @pytest.mark.parametrize(['xarray_type', 'grid_kwargs'], [(XARRAY_TYPES[0], GRID_KWARGS[0])], indirect=True)
    def test_integration_blocks(self, tmp_path, xarray_type, grid_kwargs, n_subdomains_blocks):
        """
        test the whole SubdomainCreater using block splitting
        """
        data_file = pl.Path(tmp_path, 'data.nc')
        xarray_type[0].to_netcdf(data_file)

        subs = SubdomainCreater(filename=data_file,
                                grid_kwargs=grid_kwargs,
                                target_dir=tmp_path,
                                n_subdomains=n_subdomains_blocks[0],
                                split_dims=('x', 'y'),
                                y=DIM_NAMES[0],
                                x=DIM_NAMES[1],
                                **xarray_type[1])

        # make sure missing values are set
        subs.split(use_mask=False)
        mask = xr.DataArray(np.zeros(shape=ARRAY_SHAPE, dtype=bool), dims=DIM_NAMES)

        for i_grid, slices in enumerate(n_subdomains_blocks[1], 1):
            sub_path = pl.Path(tmp_path, DEFAULT_FILE_PATTERN.format(i_grid=i_grid))
            assert sub_path.exists()
            sub = read_netcdf_file(sub_path, **xarray_type[1])
            mask[...] = False
            mask[slices] = True
            assert sub.equals(xarray_type[0].where(mask, drop=True))

    @pytest.mark.parametrize(['xarray_type', 'grid_kwargs'], [(XARRAY_TYPES[0], GRID_KWARGS[0])], indirect=True)
    def test_integration_columns(self, tmp_path, xarray_type, grid_kwargs, n_subdomains_columns):
        """
        test the whole SubdomainCreater using block splitting
        """
        data_file = pl.Path(tmp_path, 'data.nc')
        xarray_type[0].to_netcdf(data_file)

        subs = SubdomainCreater(filename=data_file,
                                grid_kwargs=grid_kwargs,
                                target_dir=tmp_path,
                                n_subdomains=n_subdomains_columns[0],
                                split_dims=('x',),
                                y=DIM_NAMES[0],
                                x=DIM_NAMES[1],
                                **xarray_type[1])

        # make sure missing values are set
        subs.split(use_mask=False)
        mask = xr.DataArray(np.zeros(shape=ARRAY_SHAPE, dtype=bool), dims=DIM_NAMES)

        for i_grid, slices in enumerate(n_subdomains_columns[1], 1):
            sub_path = pl.Path(tmp_path, DEFAULT_FILE_PATTERN.format(i_grid=i_grid))
            assert sub_path.exists()
            sub = read_netcdf_file(sub_path, **xarray_type[1])
            mask[...] = False
            mask[slices] = True
            assert sub.equals(xarray_type[0].where(mask, drop=True))

    @pytest.mark.parametrize(['xarray_type', 'grid_kwargs'], [(XARRAY_TYPES[0], GRID_KWARGS[0])], indirect=True)
    def test_integration_dxdy(self, tmp_path, xarray_type, grid_kwargs, dxdy):
        """
        test the whole SubdomainCreater using given width splitting
        """
        data_file = pl.Path(tmp_path, 'data.nc')
        xarray_type[0].to_netcdf(data_file)

        subs = SubdomainCreater(filename=data_file,
                                grid_kwargs=grid_kwargs,
                                target_dir=tmp_path,
                                dxdy=dxdy[0],
                                split_dims=('x', 'y'),
                                y=DIM_NAMES[0],
                                x=DIM_NAMES[1],
                                **xarray_type[1])

        # make sure missing values are set
        subs.split(use_mask=False)
        mask = xr.DataArray(np.zeros(shape=ARRAY_SHAPE, dtype=bool), dims=DIM_NAMES)

        for i_grid, slices in enumerate(dxdy[1], 1):
            sub_path = pl.Path(tmp_path, DEFAULT_FILE_PATTERN.format(i_grid=i_grid))
            assert sub_path.exists()
            sub = read_netcdf_file(sub_path, **xarray_type[1])
            mask[...] = False
            mask[slices] = True
            assert sub.equals(xarray_type[0].where(mask, drop=True))

    def test_check_ll(self):
        area_ll = (3, 7)
        dxdy = (2, 3)
        lls = [(3, 7), (5, 7), (-1, 7), (3, -2), (3, 10), (3, 3), (3, 12), (-2, 7), (800, 7),
               (3.0, 7.0), (5.0, 7.0), (-1.0, 7.0), (3.0, -2.0), (3.0, 10.0), (3.0, 2.4), (3.0, 12.1), (-2.5, 7.0),
               (800.199, 7.0)]
        ref_lls = [(3, 7), (3, 7), (3, 7), (3, 7), (3, 7), (3, 6), (3, 6), (2, 7), (2, 7),
                   (3.0, 7.0), (3.0, 7.0), (3.0, 7.0), (3.0, 7.0), (3.0, 7.0), (3.0, 5.4), (3.0, 6.1), (1.5, 7.0),
                   (2.199, 7.0)]
        for ref_ll, ll in zip(ref_lls, lls):
            assert np.isclose(ref_ll, check_ll(ll, dxdy, area_ll)).all()

    def test_concat(self):
        """
        - load data from multiple netcdf files
        - put Grids in a specific container
        - combine data
          - specify additional slices of other dimensions to loop on
        - dump global data
        """
        assert True
