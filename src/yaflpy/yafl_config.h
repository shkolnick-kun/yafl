/*******************************************************************************
    Copyright 2021 anonimous <shkolnick-kun@gmail.com> and contributors.

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

#ifndef YAFL_CONFIG_H
#define YAFL_CONFIG_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define YAFL_DBG(...) fprintf(stderr, __VA_ARGS__)

typedef int32_t yaflInt;

#ifndef YAFL_USE_64_BIT
#   define YAFL_USE_64_BIT (0)
#endif/*YAFL_USE_64_BIT*/

#if YAFL_USE_64_BIT
    typedef double  yaflFloat;
#   define YAFL_EPS  (1.0e-15)
#   define YAFL_SQRT sqrt
#   define YAFL_ABS  fabs
#   define YAFL_EXP  exp
#   define YAFL_LOG  log
#else/*YAFL_USE_64_BIT*/
    typedef float  yaflFloat;
#   define YAFL_EPS  (1.0e-6)
#   define YAFL_SQRT sqrtf
#   define YAFL_ABS  fabsf
#   define YAFL_EXP  expf
#   define YAFL_LOG  logf
#endif/*YAFL_USE_64_BIT*/

#ifdef __GNUC__
#   define YAFL_UNLIKELY(x) __builtin_expect((x), 0)
#else
#   define YAFL_UNLIKELY(x) (x)
#endif

#endif // YAFL_CONFIG_H

