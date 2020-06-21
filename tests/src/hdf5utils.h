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

#ifndef HDF5UTILS_H_
#define HDF5UTILS_H_

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <hdf5/serial/hdf5.h>

#define HDF5UTILS_DTYPE  H5T_FLOAT
#define HDF5UTILS_DTSIZE 8
#define HDF5UTILS_RTYPE  H5T_NATIVE_DOUBLE

typedef union {
    /*I need only up to three dims*/
    hsize_t dims[3];
    /*Named ones*/
    struct {
        hsize_t x;
        hsize_t y;
        hsize_t z;
    } dim;
} hdf5UtilsShapeSt;

typedef struct {
    hdf5UtilsShapeSt shape;
    double * data;
} hdf5UtilsMatSt;

hdf5UtilsMatSt hdf5_utils_read_array(hid_t file, const char * dsname);

void hdf5_utils_write_array(hid_t file, const char * dsname, hdf5UtilsMatSt * mat);

#endif // HDF5UTILS_H_
