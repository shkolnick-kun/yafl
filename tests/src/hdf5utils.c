/*******************************************************************************
    Copyright 2020 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
    and limitations under the License.
******************************************************************************/

#include <assert.h>
#include "hdf5utils.h"

hdf5UtilsMatSt hdf5_utils_read_array(hid_t file, const char * dsname)
{
    /* Handles */
    hid_t        space;
    hid_t        dset;
    hid_t        dcpl;
    hid_t        dtype;
    /* */
    H5D_layout_t layout;
    herr_t       status;
    int          n;
    /**/
    hdf5UtilsMatSt result = {
        .shape = {
            .dims = {0, 0, 0},
        },
        .data = NULL
    };

    dset = H5Dopen(file, dsname, H5P_DEFAULT);
    assert(dset >= 0);

    dcpl = H5Dget_create_plist(dset);
    assert(dcpl >= 0);

    layout = H5Pget_layout(dcpl);
    assert(H5D_CONTIGUOUS == layout);

    space = H5Dget_space(dset);
    assert(space >= 0);
    assert(H5Sget_simple_extent_type(space) == H5S_SIMPLE);

    n = H5Sget_simple_extent_ndims(space);
    assert(n > 0);
    assert(n <= 3);

    n = H5Sget_simple_extent_dims(space, result.shape.dims, NULL);
    assert(n > 0);

    dtype = H5Dget_type(dset);
    assert(HDF5UTILS_DTYPE == H5Tget_class(dtype));
    assert(HDF5UTILS_DTSIZE == H5Tget_size(dtype));

    n = H5Dget_storage_size(dset);
    assert(n > 0);

    result.data = (double *)malloc((size_t)n);
    assert(NULL != result.data);

    status = H5Dread(dset, HDF5UTILS_RTYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data);
    assert(status >= 0);

    status = H5Tclose(dtype);
    status = H5Sclose(space);
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    return result;
}

void hdf5_utils_write_array(hid_t file, const char * dsname, hdf5UtilsMatSt * mat)
{
    /* Handles */
    hid_t        space;
    hid_t        dset;
    hid_t        dcpl;
    hid_t        dtype;

    herr_t       status;

    int i;

    assert(mat);
    assert(mat->data);
    /*Check dims*/
    assert(mat->shape.dim.x > 0);
    if (mat->shape.dim.z > 0)
    {
        assert(mat->shape.dim.y > 0);
    }

    for (i = 0; i < 3; i++)
    {
        if (0 == mat->shape.dims[i])
        {
            break;
        }
    }
    space = H5Screate_simple(i, mat->shape.dims, NULL);
    assert(space >= 0);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    assert(dcpl >= 0);

    status = H5Pset_layout(dcpl, H5D_CONTIGUOUS);
    assert(status >= 0);

    dtype = H5Tcopy(HDF5UTILS_RTYPE);
    assert(dtype >= 0);

    dset = H5Dcreate(file, dsname, dtype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    assert(dset >= 0);

    status = H5Dwrite (dset, HDF5UTILS_RTYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat->data);
    assert(status >= 0);

    status = H5Tclose(dtype);
    status = H5Sclose(space);
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
}
