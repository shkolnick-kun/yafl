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

#ifndef YAKF_CONFIG_H
#define YAKF_CONFIG_H

#include <assert.h>
#include <math.h>
#include <stdint.h>

#define YAKF_ASSERT assert

typedef double  yakfFloat;
typedef int32_t yakfInt;

/*TODO: уточнить*/
#define YAKF_UDU_EPS  (1.0e-15)

#define YAKF_SQRT sqrt

#endif // YAKF_CONFIG_H
