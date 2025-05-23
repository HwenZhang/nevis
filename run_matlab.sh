#!/bin/bash

# 列出所有要运行的 MATLAB 脚本
# scripts=("nevis_regional_RACMO_alpha0_05_kls1e0_mu1e1_c00" \
#          "nevis_regional_RACMO_alpha0_05_kls1e0_mu1e1_c01e_1" \
#          "nevis_regional_RACMO_alpha0_05_kls1e0_mu1e1_c01e_2")

scripts=("nreg_RACMO_cg0_00_a0_2_k1_mu1e1_c1_V1e8_t300" \
         "nreg_RACMO_cg0_02_a0_2_k0_mu1e1_c1_V1e8_t300" \
         "nreg_RACMO_cg0_02_a0_2_k1_mu1e1_c1_V1e8_t300" \
         "nreg_RACMO_cg0_00_a0_05_k0_mu1e1_c1_V1e8_t300" \
         "nreg_RACMO_cg0_00_a0_05_k1_mu1e1_c1_V1e8_t300")

# 循环并并行运行每个 MATLAB 脚本
for script in "${scripts[@]}"
do
    echo "Launching ${script}"
    matlab -batch "${script}" > temp_${script}.log 2>&1

    # 检查是否执行成功
    if [ $? -ne 0 ]; then
        echo "Error running ${script}, stopping execution."
        exit 1
    fi

done

echo "All MATLAB scripts completed."