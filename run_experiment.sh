#!/bin/bash
# =============================================================================
# Efficiently runs a series of barotropic model experiments. Normal modes
# are calculated once per damping configuration, then optimal modes are
# calculated for multiple tau values.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================


# --- Configuration ---
season="JA"
echo "--- Setting season to $season ---"
sed -i "s/^season = .*/season = \"$season\"/" model_parameters.py

# Define the optimization times (tau) to loop through.
declare -a taus=(
    "2.0" "6.0" "10.0" "14.0" "18.0"
)

# Define the damping parameter sets to loop through.
# Format: "damp_flag damp_ocean damp_land"
declare -a dampings=(
    "1 20.0 10.0"
    "1 10.0 5.0"
)

# --- Functions ---
# Updates the physical damping parameters.
update_damping_params() {
    echo "Updating damping parameters: flag=$1, ocean=$2, land=$3"
    sed -i "s/^damp_flag = .*/damp_flag = $1/" model_parameters.py
    sed -i "s/^damp_ocean = .*/damp_ocean = $2/" model_parameters.py
    sed -i "s/^damp_land = .*/damp_land = $3/" model_parameters.py
}

# Updates only the tau parameter.
update_tau() {
    echo "Updating tau = $1"
    sed -i "s/^tau = .*/tau = $1/" model_parameters.py
}

# Get normal modes.
run_normal_mode_analysis() {
    echo "Running normal mode analysis..."
    python 00_damping_upper.py
    python 01_LHS_Matrix.py
    python 02_RHS_Matrix.py
    python 03_normal_mode.py
    python plot_normal_mode.py
}

# Get optimal modes for the current tau.
run_optimal_mode_analysis() {
    echo "Running optimal mode analysis for tau..."
    python 04_optimal_mode.py
    python plot_optimal_mode.py
}

# --- Main Execution Loop ---
for damp_params in "${dampings[@]}"; do
    read -r d_flag d_ocean d_land <<< "$damp_params"
    
    echo -e "\n=============================================================================="
    echo "Processing Damping Config: flag=$d_flag, ocean=$d_ocean, land=$d_land"
    
    update_damping_params "$d_flag" "$d_ocean" "$d_land"
    run_normal_mode_analysis
    
    for tau_val in "${taus[@]}"; do
        echo "------------------------------------------------------------------------------"
        echo "Calculating optimal modes for tau = $tau_val"
        update_tau "$tau_val"
        run_optimal_mode_analysis
    done
done

echo -e "\n--- All experiments completed successfully. ---"
