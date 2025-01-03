import streamlit as st
import pandas as pd
import numpy as np
import altair as alt

#############################
#  1) Core Calculation
#############################
def compute_prmt_distribution(M, Kd_PM, S, Kd_PS, I, Kd_PI, Kd_PMI, Kd_PSI):
    """
    Computes fractional occupancy of a PRMT among 6 microstates:
      1) free (P_free)
      2) bound to MTA (P_M)
      3) bound to SAM (P_S)
      4) bound to inhibitor (P_I)
      5) bound to MTA + inhibitor (P_MI)
      6) bound to SAM + inhibitor (P_SI)

    All Kd and concentration values must be in consistent units (e.g. µM).

    Returns a dict with fractions for each microstate.
    """
    alpha_M  = M / Kd_PM   # dimensionless ratio for MTA binding
    alpha_S  = S / Kd_PS   # ratio for SAM binding
    alpha_Ip = I / Kd_PI   # ratio for inhibitor binding free protein
    alpha_IM = I / Kd_PMI  # ratio for inhibitor binding MTA-bound protein
    alpha_IS = I / Kd_PSI  # ratio for inhibitor binding SAM-bound protein

    # Microstate weights
    W_free = 1.0
    W_PM   = alpha_M
    W_PS   = alpha_S
    W_PI   = alpha_Ip
    W_PMI  = alpha_M * alpha_IM
    W_PSI  = alpha_S * alpha_IS

    Z = W_free + W_PM + W_PS + W_PI + W_PMI + W_PSI

    return {
        "P_free": W_free / Z,
        "P_M":    W_PM   / Z,
        "P_S":    W_PS   / Z,
        "P_I":    W_PI   / Z,
        "P_MI":   W_PMI  / Z,
        "P_SI":   W_PSI  / Z,
    }


#############################
#  2) Main Streamlit App
#############################
def main():
    st.title("Interactive PRMT Inhibition App (Log-scale x-axis)")

    st.markdown(
        """
        This app computes the *functional fraction* of PRMT (i.e., \\(P_{\\mathrm{free}} + P_{S}\\)) 
        under two cell conditions (**Healthy** vs. **MTAP-deleted**), for both **PRMT5** and **PRMT1**. 
        It uses **Altair** to display the x-axis in *log scale*.  
        
        **Instructions**:  
        1. Adjust the parameters in the sidebar (MTA, SAM, Kd values, etc.).  
        2. A range of inhibitor concentrations (in nM) is swept on a *log scale*.  
        3. See how the “functional fraction” curves shift for PRMT5 and PRMT1!
        """
    )

    # --- Sidebar Inputs ---
    with st.sidebar:
        st.header("Cell Conditions")
        MTA_healthy = st.number_input("MTA (Healthy) [µM]", min_value=0.0, value=1.0, step=0.1)
        MTA_deleted = st.number_input("MTA (MTAP-deleted) [µM]", min_value=0.0, value=100.0, step=1.0)
        SAM_conc    = st.number_input("SAM [µM]", min_value=0.0, value=50.0, step=1.0)
        
        st.header("PRMT5 Kd (µM)")
        Kd_PM_5 = st.number_input("Kd(PRMT5–MTA)", min_value=0.0001, value=5.0, step=0.1)
        Kd_PS_5 = st.number_input("Kd(PRMT5–SAM)", min_value=0.0001, value=10.0, step=0.1)
        
        st.header("PRMT1 Kd (µM)")
        Kd_PM_1 = st.number_input("Kd(PRMT1–MTA)", min_value=0.0001, value=50.0, step=1.0)
        Kd_PS_1 = st.number_input("Kd(PRMT1–SAM)", min_value=0.0001, value=10.0, step=0.1)

        st.header("Inhibitor Kd (µM)")
        Kd_PI   = st.number_input("Kd(Inhib–free)", min_value=0.000001, value=0.004, step=0.001, format="%.6f")
        Kd_PMI  = st.number_input("Kd(Inhib–MTA-bound)", min_value=0.000001, value=0.004, step=0.001, format="%.6f")
        Kd_PSI  = st.number_input("Kd(Inhib–SAM-bound)", min_value=0.000001, value=0.004, step=0.001, format="%.6f")

        st.header("Inhibitor Range (nM)")
        min_inhib_nM = st.number_input("Min [Inhibitor] (nM)", min_value=0.001, value=1.0, step=0.5, format="%.3f")
        max_inhib_nM = st.number_input("Max [Inhibitor] (nM)", min_value=1.0, value=1e4, step=1000.0, format="%.1f")
        num_points   = st.slider("Number of log points", min_value=10, max_value=100, value=50)

    # --- Generate inhibitor array in nM (log space) ---
    I_values_nM = np.logspace(
        np.log10(min_inhib_nM),
        np.log10(max_inhib_nM),
        num_points
    )

    # --- Prepare PRMT5 data (wide format) ---
    prmt5_data = {
        "Inhibitor (nM)": I_values_nM,
        "Healthy": [],
        "MTAP-deleted": []
    }

    # --- Prepare PRMT1 data (wide format) ---
    prmt1_data = {
        "Inhibitor (nM)": I_values_nM,
        "Healthy": [],
        "MTAP-deleted": []
    }

    # --- Calculate functional fraction for each row ---
    for InM in I_values_nM:
        I_uM = InM * 1e-3  # Convert nM -> µM

        # PRMT5 - healthy
        fracs_5_h = compute_prmt_distribution(
            M=MTA_healthy, Kd_PM=Kd_PM_5,
            S=SAM_conc,    Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_5_h = fracs_5_h["P_free"] + fracs_5_h["P_S"]
        prmt5_data["Healthy"].append(func_5_h)

        # PRMT5 - MTAP-del
        fracs_5_d = compute_prmt_distribution(
            M=MTA_deleted, Kd_PM=Kd_PM_5,
            S=SAM_conc,    Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_5_d = fracs_5_d["P_free"] + fracs_5_d["P_S"]
        prmt5_data["MTAP-deleted"].append(func_5_d)

        # PRMT1 - healthy
        fracs_1_h = compute_prmt_distribution(
            M=MTA_healthy, Kd_PM=Kd_PM_1,
            S=SAM_conc,    Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_1_h = fracs_1_h["P_free"] + fracs_1_h["P_S"]
        prmt1_data["Healthy"].append(func_1_h)

        # PRMT1 - MTAP-del
        fracs_1_d = compute_prmt_distribution(
            M=MTA_deleted, Kd_PM=Kd_PM_1,
            S=SAM_conc,    Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_1_d = fracs_1_d["P_free"] + fracs_1_d["P_S"]
        prmt1_data["MTAP-deleted"].append(func_1_d)

    # --- Convert to DataFrame (wide format) ---
    df_prmt5 = pd.DataFrame(prmt5_data)
    df_prmt1 = pd.DataFrame(prmt1_data)

    ###########################
    #  3) Altair Log-Scale Plot
    ###########################
    st.subheader("PRMT5 Results (Log-scale x-axis)")

    # Melt wide -> long
    df5_long = df_prmt5.melt(
        id_vars="Inhibitor (nM)", 
        var_name="Condition",
        value_name="Fraction"
    )

    chart5 = (
        alt.Chart(df5_long)
        .mark_line()
        .encode(
            x=alt.X(
                "Inhibitor (nM):Q",
                title="[Inhibitor] (nM)",
                scale=alt.Scale(type="log")
            ),
            y=alt.Y(
                "Fraction:Q",
                title="Fraction (P_free + P_S)",
                scale=alt.Scale(domain=[0, 1])  # fraction range
            ),
            color=alt.Color("Condition:N", legend=alt.Legend(title="Cell Type"))
        )
        .properties(width="container", height=400)
    )
    st.altair_chart(chart5, use_container_width=True)

    st.subheader("PRMT1 Results (Log-scale x-axis)")

    # Melt wide -> long
    df1_long = df_prmt1.melt(
        id_vars="Inhibitor (nM)",
        var_name="Condition",
        value_name="Fraction"
    )

    chart1 = (
        alt.Chart(df1_long)
        .mark_line()
        .encode(
            x=alt.X(
                "Inhibitor (nM):Q",
                title="[Inhibitor] (nM)",
                scale=alt.Scale(type="log")
            ),
            y=alt.Y(
                "Fraction:Q",
                title="Fraction (P_free + P_S)",
                scale=alt.Scale(domain=[0, 1])
            ),
            color=alt.Color("Condition:N", legend=alt.Legend(title="Cell Type"))
        )
        .properties(width="container", height=400)
    )
    st.altair_chart(chart1, use_container_width=True)

    st.markdown("---")
    st.markdown("**Note**: The x-axis is in nM (log scale). You can adjust the range of inhibitor concentrations and other parameters in the sidebar, and the charts will update automatically.")

if __name__ == "__main__":
    main()

