#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright 2022 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
"""
import time
import matplotlib.pyplot as plt
import numpy as np

from ab_tests import *
from case1    import *

from yaflpy import Joseph
from yaflpy import Bierman
from yaflpy import AdaptiveBierman
from yaflpy import AdaptiveJoseph
from yaflpy import RobustBierman
from yaflpy import RobustJoseph
from yaflpy import AdaptiveRobustBierman
from yaflpy import AdaptiveRobustJoseph

from yaflpy import Unscented
from yaflpy import UnscentedAdaptive
from yaflpy import UnscentedAdaptiveBierman
from yaflpy import UnscentedAdaptiveRobustBierman
from yaflpy import UnscentedBierman
from yaflpy import UnscentedRobustBierman
#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clear, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
a = case_ekf(Joseph, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/j_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(Bierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/b_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(AdaptiveJoseph, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/aj_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(AdaptiveBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/ab_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(RobustJoseph, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/rj_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(RobustBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/rb_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(AdaptiveRobustJoseph, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/arj_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ekf(AdaptiveRobustBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/arb_ekf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(Unscented, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/ukf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(UnscentedAdaptive, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/aukf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(UnscentedBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/b_ukf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(UnscentedAdaptiveBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/ab_ukf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(UnscentedRobustBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/rb_ukf_case1_64bit.h5')

#------------------------------------------------------------------------------
a = case_ukf(UnscentedAdaptiveRobustBierman, STD)
yafl_record_test_data(a, clear, noisy, t, '../../data/arb_ukf_case1_64bit.h5')
