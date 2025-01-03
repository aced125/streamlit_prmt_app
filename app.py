import streamlit as st
import numpy as np
import pandas as pd

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
#  2) Streamlit UI
#############################
def main():
    st.title("Interactive PRMT Inhibition App")

    st.markdown(
        """
        This app computes the *functional fraction* of PRMT (i.e., \\(P_{\\mathrm{free}} + P_{S}\\)) 
        under two cell conditions (Healthy vs. MTAP-deleted), for **both PRMT5 and PRMT1**. 
        Adjust the parameters in the sidebar to see how the curves change!
        """
    )

    with st.sidebar:
        st.header("Cell Conditions")
        MTA_healthy = st.number_input("MTA (Healthy) [µM]", min_value=0.0, value=1.0, step=0.1)
        MTA_deleted = st.number_input("MTA (MTAP-deleted) [µM]", min_value=0.0, value=100.0, step=1.0)
        SAM_conc    = st.number_input("SAM [µM]", min_value=0.0, value=50.0, step=1.0)
        
        st.header("PRMT5 Kd (µM)")
        Kd_PM_5 = st.number_input("Kd(PRMT5–MTA) [µM]", min_value=0.0001, value=5.0, step=0.1)
        Kd_PS_5 = st.number_input("Kd(PRMT5–SAM) [µM]", min_value=0.0001, value=10.0, step=0.1)
        
        st.header("PRMT1 Kd (µM)")
        Kd_PM_1 = st.number_input("Kd(PRMT1–MTA) [µM]", min_value=0.0001, value=50.0, step=1.0)
        Kd_PS_1 = st.number_input("Kd(PRMT1–SAM) [µM]", min_value=0.0001, value=10.0, step=0.1)

        st.header("Inhibitor Kd (µM)")
        Kd_PI   = st.number_input("Kd(Inhib–free) [µM]", min_value=0.000001, value=0.004, step=0.001, format="%.6f")
        Kd_PMI  = st.number_input("Kd(Inhib–MTA-bound) [µM]", min_value=0.000001, value=0.004, step=0.001, format="%.6f")
        Kd_PSI  = st.number_input("Kd(Inhib–SAM-bound) [µM]", min_value=0.000001, value=0.004, step=0.001, format="%.6f")

        st.header("Inhibitor Range (nM)")
        min_inhib_nM = st.number_input("Min [Inhibitor] (nM)", min_value=0.001, value=1.0, step=0.5, format="%.3f")
        max_inhib_nM = st.number_input("Max [Inhibitor] (nM)", min_value=1.0, value=1e4, step=1000.0, format="%.1f")
        num_points   = st.slider("Number of points in log range", min_value=10, max_value=100, value=50)
    
    # Create array of inhibitor concentrations in nM
    I_values_nM = np.logspace(np.log10(min_inhib_nM), np.log10(max_inhib_nM), num_points)

    # We'll build data frames for PRMT5 and PRMT1 to show in line charts.
    # Each data frame has columns:
    #   "Inhibitor (nM)", "Healthy", "MTAP-deleted"
    # so we can do st.line_chart(...) specifying y=["Healthy", "MTAP-deleted"] if we want.
    prmt5_data = {
        "Inhibitor (nM)": I_values_nM,
        "Healthy": [],
        "MTAP-deleted": []
    }
    prmt1_data = {
        "Inhibitor (nM)": I_values_nM,
        "Healthy": [],
        "MTAP-deleted": []
    }

    for InM in I_values_nM:
        I_uM = InM * 1e-3  # convert nM -> µM

        # PRMT5 - healthy
        fracs_5_healthy = compute_prmt_distribution(
            M=MTA_healthy, Kd_PM=Kd_PM_5,
            S=SAM_conc,    Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_5_healthy = fracs_5_healthy["P_free"] + fracs_5_healthy["P_S"]
        prmt5_data["Healthy"].append(func_5_healthy)

        # PRMT5 - MTAP-deleted
        fracs_5_deleted = compute_prmt_distribution(
            M=MTA_deleted, Kd_PM=Kd_PM_5,
            S=SAM_conc,    Kd_PS=Kd_PS_5,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_5_deleted = fracs_5_deleted["P_free"] + fracs_5_deleted["P_S"]
        prmt5_data["MTAP-deleted"].append(func_5_deleted)

        # PRMT1 - healthy
        fracs_1_healthy = compute_prmt_distribution(
            M=MTA_healthy, Kd_PM=Kd_PM_1,
            S=SAM_conc,    Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_1_healthy = fracs_1_healthy["P_free"] + fracs_1_healthy["P_S"]
        prmt1_data["Healthy"].append(func_1_healthy)

        # PRMT1 - MTAP-deleted
        fracs_1_deleted = compute_prmt_distribution(
            M=MTA_deleted, Kd_PM=Kd_PM_1,
            S=SAM_conc,    Kd_PS=Kd_PS_1,
            I=I_uM,
            Kd_PI=Kd_PI,   Kd_PMI=Kd_PMI, Kd_PSI=Kd_PSI
        )
        func_1_deleted = fracs_1_deleted["P_free"] + fracs_1_deleted["P_S"]
        prmt1_data["MTAP-deleted"].append(func_1_deleted)

    # Convert to DataFrame
    df_prmt5 = pd.DataFrame(prmt5_data)
    df_prmt1 = pd.DataFrame(prmt1_data)

    st.subheader("PRMT5 Results")
    st.write("Fraction of PRMT5 that is free + SAM-bound vs. [Inhibitor] (nM)")
    st.line_chart(
        data=df_prmt5,
        x="Inhibitor (nM)",
        y=["Healthy", "MTAP-deleted"],
        use_container_width=True
    )

    st.subheader("PRMT1 Results")
    st.write("Fraction of PRMT1 that is free + SAM-bound vs. [Inhibitor] (nM)")
    st.line_chart(
        data=df_prmt1,
        x="Inhibitor (nM)",
        y=["Healthy", "MTAP-deleted"],
        use_container_width=True
    )

    st.markdown("---")
    st.markdown("**Tip**: Try changing the parameters in the sidebar (Kd values, cell conditions, etc.) to see how the curves move in real time!")

if __name__ == "__main__":
    main()

