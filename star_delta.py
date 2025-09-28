import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

st.set_page_config(page_title="Star–Delta Simulation", layout="wide")
st.title("⭐🔺 Star–Delta Connection: Interactive Virtual Lab")

# -------------------- Sidebar Controls --------------------
with st.sidebar:
    st.header("Inputs")
    VL_rms = st.number_input("Supply Line Voltage VL (rms, V)", 100.0, 1000.0, 415.0, step=5.0)
    f = st.slider("Frequency (Hz)", 40, 70, 50)
    connection = st.radio("Connection", ("Star (Y)", "Delta (Δ)"))
    load_type = st.selectbox("Load Type (per phase)", ["R", "R–L", "R–C"])
    R = st.slider("Resistance R (Ω)", 1, 200, 30)
    L_mH = st.slider("Inductance L (mH)", 1, 500, 100) if load_type == "R–L" else 0
    C_uF = st.slider("Capacitance C (µF)", 1, 500, 50) if load_type == "R–C" else 0
    t_end_ms = st.slider("Waveform window (ms)", 5, 40, 20)
    snapshot_deg = st.slider("Phasor snapshot angle (°)", 0, 359, 0)

# -------------------- Derived Quantities --------------------
w = 2*np.pi*f
L = L_mH/1000.0      # H
C = C_uF*1e-6        # F
Xl = w*L             # Ω
Xc = (1/(w*C)) if (load_type == "R–C" and C > 0) else 0.0

# Per-phase voltage (rms)
if connection == "Star (Y)":
    Vph_rms = VL_rms/np.sqrt(3)
else:  # Delta
    Vph_rms = VL_rms

# Per-phase impedance
if load_type == "R":
    Zph = R + 0j
elif load_type == "R–L":
    Zph = complex(R, Xl)
else:                      # R–C
    Zph = complex(R, -Xc) if Xc != 0 else complex(R, 0)

# Per-phase current (phasor, rms)
Iph = (Vph_rms / Zph) if Zph != 0 else complex(np.inf, 0)
Iph_mag = np.abs(Iph)
Iph_ang = np.degrees(np.angle(Iph))
phi_Z = np.angle(Zph) if Zph != 0 else 0.0
pf = np.cos(phi_Z) if Zph != 0 else 1.0
pf_type = "lagging" if (load_type == "R–L") else ("leading" if (load_type == "R–C") else "unity")

# Line current magnitude (rms)
IL_mag = Iph_mag if connection == "Star (Y)" else np.sqrt(3)*Iph_mag

# Power (three-phase)
P_phase = Vph_rms * Iph_mag * pf
P_total = 3 * P_phase
Q_phase = Vph_rms * Iph_mag * np.sin(phi_Z)
Q_total = 3 * Q_phase
S_total = np.sqrt(P_total**2 + Q_total**2)

# -------------------- KPIs --------------------
c1, c2, c3, c4 = st.columns(4)
c1.metric("Phase Voltage Vph (rms)", f"{Vph_rms:.1f} V")
c2.metric("Line Current IL (rms)", f"{IL_mag:.2f} A")
c3.metric("Power Factor", f"{abs(pf):.3f} ({pf_type})")
c4.metric("Total Real Power P", f"{P_total/1000:.2f} kW")

# -------------------- Waveforms --------------------
t = np.linspace(0, t_end_ms/1000.0, 2000)

# Phase-to-neutral voltages (instantaneous)
Va = np.sqrt(2)*Vph_rms*np.sin(w*t)
Vb = np.sqrt(2)*Vph_rms*np.sin(w*t - 2*np.pi/3)
Vc = np.sqrt(2)*Vph_rms*np.sin(w*t + 2*np.pi/3)

# Line-to-line voltages (instantaneous): Vab = Va - Vb etc.
Vab = Va - Vb
Vbc = Vb - Vc
Vca = Vc - Va

# Instantaneous phase currents (shift by -angle(Z))
phi_I = -phi_Z
Ia = np.sqrt(2)*Iph_mag*np.sin(w*t + phi_I)
Ib = np.sqrt(2)*Iph_mag*np.sin(w*t - 2*np.pi/3 + phi_I)
Ic = np.sqrt(2)*Iph_mag*np.sin(w*t + 2*np.pi/3 + phi_I)

tab1, tab2, tab3, tab4 = st.tabs(["📈 Waveforms", "🧭 Phasors", "🧪 Quiz (3 Qs)", "📄 Lab Manual & Export"])

with tab1:
    colA, colB = st.columns(2)

    # Phase voltages
    fig_vp, ax_vp = plt.subplots()
    ax_vp.plot(t*1000, Va, label="Va")
    ax_vp.plot(t*1000, Vb, label="Vb")
    ax_vp.plot(t*1000, Vc, label="Vc")
    ax_vp.set_xlabel("Time (ms)")
    ax_vp.set_ylabel("Voltage (V)")
    ax_vp.set_title("Phase Voltages (instantaneous)")
    ax_vp.legend()
    colA.pyplot(fig_vp)

    # Line-to-line voltages
    fig_vl, ax_vl = plt.subplots()
    ax_vl.plot(t*1000, Vab, label="Vab")
    ax_vl.plot(t*1000, Vbc, label="Vbc")
    ax_vl.plot(t*1000, Vca, label="Vca")
    ax_vl.set_xlabel("Time (ms)")
    ax_vl.set_ylabel("Voltage (V)")
    ax_vl.set_title("Line-to-Line Voltages (instantaneous)")
    ax_vl.legend()
    colB.pyplot(fig_vl)

    # Phase currents
    fig_i, ax_i = plt.subplots()
    ax_i.plot(t*1000, Ia, label="Ia")
    ax_i.plot(t*1000, Ib, label="Ib")
    ax_i.plot(t*1000, Ic, label="Ic")
    ax_i.set_xlabel("Time (ms)")
    ax_i.set_ylabel("Current (A)")
    ax_i.set_title("Phase Currents (instantaneous)")
    ax_i.legend()
    st.pyplot(fig_i)

with tab2:
    theta = np.radians(snapshot_deg)
    # Voltage phasors (rms)
    Vphasors = [
        (Vph_rms*np.cos(theta), Vph_rms*np.sin(theta), "Va"),
        (Vph_rms*np.cos(theta-2*np.pi/3), Vph_rms*np.sin(theta-2*np.pi/3), "Vb"),
        (Vph_rms*np.cos(theta+2*np.pi/3), Vph_rms*np.sin(theta+2*np.pi/3), "Vc"),
    ]
    # Current phasors (rms): shifted by phi_I
    Iphasors = [
        (Iph_mag*np.cos(theta+phi_I), Iph_mag*np.sin(theta+phi_I), "Ia"),
        (Iph_mag*np.cos(theta-2*np.pi/3+phi_I), Iph_mag*np.sin(theta-2*np.pi/3+phi_I), "Ib"),
        (Iph_mag*np.cos(theta+2*np.pi/3+phi_I), Iph_mag*np.sin(theta+2*np.pi/3+phi_I), "Ic"),
    ]

    fig_p, ax_p = plt.subplots()
    for x, y, lbl in Vphasors:
        ax_p.arrow(0, 0, x, y, head_width=0.05*Vph_rms, length_includes_head=True)
        ax_p.text(x*1.05, y*1.05, lbl)
    for x, y, lbl in Iphasors:
        ax_p.arrow(0, 0, x, y, head_width=0.05*Iph_mag, length_includes_head=True)
        ax_p.text(x*1.05, y*1.05, lbl)

    vmax = max(Vph_rms, Iph_mag)
    ax_p.set_xlim(-1.6*vmax, 1.6*vmax)
    ax_p.set_ylim(-1.6*vmax, 1.6*vmax)
    ax_p.set_aspect("equal", "box")
    ax_p.grid(True, alpha=0.4)
    ax_p.set_title("Phasor Diagram (rms) — rotate with snapshot angle")
    ax_p.set_xlabel("Real axis")
    ax_p.set_ylabel("Imag axis")
    st.pyplot(fig_p)

with tab3:
    st.subheader("Quick Quiz (auto-scored)")
    st.caption("Answer all three and click **Check Score**.")
    q1 = st.radio("1) In Star (Y):", [
        "V_L = V_ph and I_L = √3·I_ph",
        "V_L = √3·V_ph and I_L = I_ph",
        "V_L = √3·V_ph and I_L = √3·I_ph"
    ], index=None)
    q2 = st.radio("2) In Delta (Δ):", [
        "V_L = V_ph and I_L = √3·I_ph",
        "V_L = √3·V_ph and I_L = I_ph",
        "V_L = V_ph and I_L = I_ph"
    ], index=None)
    q3 = st.radio("3) With R–L load, power factor is:", [
        "leading", "lagging", "unity"
    ], index=None)

    if st.button("Check Score"):
        score = 0
        if q1 == "V_L = √3·V_ph and I_L = I_ph": score += 1
        if q2 == "V_L = V_ph and I_L = √3·I_ph": score += 1
        if q3 == "lagging": score += 1
        st.success(f"Your score: {score}/3")
        if score < 3:
            st.info("Tip: switch between Star and Delta in the app and observe PF with R–L / R–C.")

with tab4:
    st.subheader("Lab Manual (auto-filled, editable)")
    principle = ("In a three-phase system, loads may be connected in Star (Y) or Delta (Δ). "
                 "In Star, line voltage is √3 times phase voltage while line current equals phase current. "
                 "In Delta, line voltage equals phase voltage while line current is √3 times phase current. "
                 "This lab visualizes voltage/current relations, phasors, and power for R, R–L, and R–C loads.")
    procedure = (
        "1) Set supply VL and frequency. Select Star or Delta.\n"
        "2) Choose load type (R/R–L/R–C) and set component values.\n"
        "3) Observe KPIs (Vph, IL, PF, P, Q, S). Rotate phasor snapshot to relate waveforms and vectors.\n"
        "4) Compare phase vs line-to-line voltages; note current phase shift with R–L (lag) and R–C (lead).\n"
        "5) Export CSV and answer the 3-question quiz.\n"
    )

    st.markdown("**Principle (≈50 words)**")
    principle_edit = st.text_area("", principle, height=80)
    st.markdown("**Instructions (≈100 words)**")
    procedure_edit = st.text_area("", procedure, height=140)

    # Show computed summary for the report
    st.markdown("**Auto-computed summary (copy to your report):**")
    st.write(f"- Connection: **{connection}** | Load: **{load_type}** with R={R} Ω, L={L_mH} mH, C={C_uF} µF")
    st.write(f"- Vph(rms)={Vph_rms:.1f} V, I_line(rms)={IL_mag:.2f} A, PF={abs(pf):.3f} ({pf_type})")
    st.write(f"- Total power: P={P_total/1000:.3f} kW, Q={Q_total/1000:.3f} kVAr, S={S_total/1000:.3f} kVA")

    # Download buttons
    N = len(t)
    data = np.column_stack([t, Va, Vb, Vc, Vab, Vbc, Vca, Ia, Ib, Ic])
    csv = "t,Va,Vb,Vc,Vab,Vbc,Vca,Ia,Ib,Ic\n" + "\n".join(
        [",".join(f"{v:.6f}" for v in row) for row in data]
    )
    st.download_button("Download waveforms (CSV)", data=csv,
                       file_name="star_delta_waveforms.csv", mime="text/csv")

    manual_text = StringIO()
    manual_text.write("# Star–Delta Virtual Lab Manual\n\n")
    manual_text.write("## Principle\n")
    manual_text.write(principle_edit + "\n\n")
    manual_text.write("## Apparatus/Software\n")
    manual_text.write("- Streamlit app (this simulation)\n- Voltmeter & Ammeter (virtual)\n- Three-phase supply (virtual)\n\n")
    manual_text.write("## Procedure\n")
    manual_text.write(procedure_edit + "\n")
    manual_text.write("## Observation Template\n")
    manual_text.write("| Conn | VL (V) | Vph (V) | IL (A) | Iph (A) | PF | P(kW) | Q(kVAr) |\n")
    manual_text.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
    manual_text.write("|\n")
    manual_text.write("## Result\nVerified star and delta relationships and observed PF behavior with R/R–L/R–C loads.\n")
    st.download_button("Download Lab Manual (text)", data=manual_text.getvalue(),
                       file_name="Star_Delta_Lab_Manual.txt", mime="text/plain")

st.caption("Reference relations: Star → V_L=√3·V_ph, I_L=I_ph | Delta → V_L=V_ph, I_L=√3·I_ph. PF lag for R–L, lead for R–C.")
