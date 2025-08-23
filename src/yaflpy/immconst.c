#include <yafl.h>

#define YAFL_CB(func)  const yaflKalmanUpdateCBP imm_##func = (yaflKalmanUpdateCBP)func
#define YAFL_CB2(func) const yaflKalmanUpdateCBP2 imm_##func = (yaflKalmanUpdateCBP2)func

YAFL_CB(yafl_ekf_base_predict);
YAFL_CB2(yafl_ekf_bierman_update);
YAFL_CB2(yafl_ekf_joseph_update);
YAFL_CB2(yafl_ekf_adaptive_bierman_update);
YAFL_CB2(yafl_ekf_adaptive_joseph_update);
YAFL_CB2(yafl_ekf_robust_bierman_update);
YAFL_CB2(yafl_ekf_robust_joseph_update);
YAFL_CB2(yafl_ekf_adaptive_robust_bierman_update);
YAFL_CB2(yafl_ekf_adaptive_robust_joseph_update);

YAFL_CB(yafl_ukf_base_predict);
YAFL_CB2(yafl_ukf_bierman_update);
YAFL_CB2(yafl_ukf_adaptive_bierman_update);
YAFL_CB2(yafl_ukf_robust_bierman_update);
YAFL_CB2(yafl_ukf_adaptive_robust_bierman_update);

YAFL_CB2(yafl_ukf_update);
YAFL_CB2(yafl_ukf_adaptive_update);
