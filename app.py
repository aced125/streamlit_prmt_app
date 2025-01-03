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

    All Kd and concentration values must be in consistent units (e.g., µM).

    Returns a dict with fractions for each microstate.
    """
    alpha_M  = M / Kd_PM   # dimensionless ratio for MTA binding
    alpha_S  = S / Kd_PS   # ratio for SAM binding
    alpha_Ip = I / Kd_PI   # ratio for inhibitor binding (protein-free)
    alpha_IM = I / Kd_PMI  # ratio for inhibitor binding (MTA-bound)
    alpha_IS = I / Kd_PSI  # ratio for inhibitor binding (SAM-bound)

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
    st.title("PRMT Inhibition App")

    st.markdown(
        """
        This app computes the *functional fraction* of PRMT (free PRMT + PRMT-SAM) 
        in **Healthy** vs. **MTAP-deleted** cells at varying inhibitor concentrations.
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

        st.header("Inhibitor Selectivity Setup")

        # 1) Kd(Inhib–MTA) in nM (number_input)
        kd_inhib_mta_nM = st.number_input(
            "Kd(Inhib–PRMT–MTA) [nM]",
            min_value=0.0001,
            value=10.0,
            step=0.5,
            help="Kd of Inhibitor to PRMT–MTA complex, in nM."
        )

        # 2) Selectivity factor (radio with 3 options)
        selectivity_factor = st.radio(
            "Selectivity factor (fold difference)",
            options=[1.0, 10.0, 100.0],
            index=1,  # default to 10
            format_func=lambda x: f"{int(x)}-fold",
            help="Ratio of Kd(Inhib–free) to Kd(Inhib–MTA). Also used for Kd(Inhib–SAM)."
        )

        st.write(f"**Example**: If Kd(Inhib–MTA) = {kd_inhib_mta_nM:.2f} nM, selectivity = {selectivity_factor}×, then:")
        st.write(f"- Kd(Inhib–free) = {kd_inhib_mta_nM*selectivity_factor:.2f} nM")
        st.write(f"- Kd(Inhib–SAM)  = {kd_inhib_mta_nM*selectivity_factor:.2f} nM")

        st.header("Inhibitor Range (nM)")
        min_inhib_nM = st.number_input("Min [Inhibitor] (nM)", min_value=0.0001, value=0.01, step=0.5, format="%.3f")
        max_inhib_nM = st.number_input("Max [Inhibitor] (nM)", min_value=1.0, value=1e4, step=1000.0, format="%.1f")
        num_points   = st.slider("Number of log points", min_value=10, max_value=100, value=50)

    # --- Convert user Kds (nM) to µM for calculations ---
    kd_inhib_mta_uM = kd_inhib_mta_nM * 1e-3  # nM -> µM
    kd_inhib_free_uM = (kd_inhib_mta_nM * selectivity_factor) * 1e-3
    kd_inhib_sam_uM  = (kd_inhib_mta_nM * selectivity_factor) * 1e-3

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
            M=MTA_healthy, 
            Kd_PM=Kd_PM_5,
            S=SAM_conc,
            Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=kd_inhib_free_uM,
            Kd_PMI=kd_inhib_mta_uM,
            Kd_PSI=kd_inhib_sam_uM
        )
        func_5_h = fracs_5_h["P_free"] + fracs_5_h["P_S"]
        prmt5_data["Healthy"].append(func_5_h)

        # PRMT5 - MTAP-del
        fracs_5_d = compute_prmt_distribution(
            M=MTA_deleted, 
            Kd_PM=Kd_PM_5,
            S=SAM_conc,    
            Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=kd_inhib_free_uM,
            Kd_PMI=kd_inhib_mta_uM,
            Kd_PSI=kd_inhib_sam_uM
        )
        func_5_d = fracs_5_d["P_free"] + fracs_5_d["P_S"]
        prmt5_data["MTAP-deleted"].append(func_5_d)

        # PRMT1 - healthy
        fracs_1_h = compute_prmt_distribution(
            M=MTA_healthy, 
            Kd_PM=Kd_PM_1,
            S=SAM_conc,
            Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=kd_inhib_free_uM,
            Kd_PMI=kd_inhib_mta_uM,
            Kd_PSI=kd_inhib_sam_uM
        )
        func_1_h = fracs_1_h["P_free"] + fracs_1_h["P_S"]
        prmt1_data["Healthy"].append(func_1_h)

        # PRMT1 - MTAP-del
        fracs_1_d = compute_prmt_distribution(
            M=MTA_deleted,
            Kd_PM=Kd_PM_1,
            S=SAM_conc,
            Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=kd_inhib_free_uM,
            Kd_PMI=kd_inhib_mta_uM,
            Kd_PSI=kd_inhib_sam_uM
        )
        func_1_d = fracs_1_d["P_free"] + fracs_1_d["P_S"]
        prmt1_data["MTAP-deleted"].append(func_1_d)

    # --- Convert to DataFrame (wide format) ---
    df_prmt5 = pd.DataFrame(prmt5_data)
    df_prmt1 = pd.DataFrame(prmt1_data)

    ###########################
    #  3) Altair Log-Scale Plot
    ###########################
    st.subheader(f"PRMT5i ({int(round(selectivity_factor))}x selective) Effect on Functional PRMT5")

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
                title="Functional PRMT5 Remaining",
                scale=alt.Scale(domain=[0, 1])  # fraction range
            ),
            color=alt.Color("Condition:N", legend=alt.Legend(title="Cell Type"))
        )
        .properties(width="container", height=400)
    )
    st.altair_chart(chart5, use_container_width=True)

    st.subheader(f"PRMT1i ({int(round(selectivity_factor))}x selective) Effect on Functional PRMT1")

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
                title="Functional PRMT1 Remaining",
                scale=alt.Scale(domain=[0, 1])
            ),
            color=alt.Color("Condition:N", legend=alt.Legend(title="Cell Type"))
        )
        .properties(width="container", height=400)
    )
    st.altair_chart(chart1, use_container_width=True)

    st.markdown("---")
    st.markdown("""
    **Notes**:
    - The slider controls *Kd(Inhib–MTA)* (in nM) and the *Selectivity factor*. 
      - If selectivity = 1, there's no difference in Kd among free, MTA-bound, or SAM-bound.
      - If selectivity = 10, the inhibitor is 10x tighter for MTA-bound vs. free or SAM-bound.
    - X-axis is log-scale in nM.
    - MTA, SAM, and the PRMT Kds are given in µM.
    """)

if __name__ == "__main__":
    main()
