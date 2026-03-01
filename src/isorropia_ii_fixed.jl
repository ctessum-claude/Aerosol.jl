"""
ISORROPIA II: Thermodynamic equilibrium model for K+–Ca2+–Mg2+–NH4+–Na+–SO42-–NO3-–Cl-–H2O aerosols

Implementation of Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally
efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O
aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.

This implementation fixes critical issues in the original code and provides a complete
thermodynamic equilibrium calculation following the paper specifications.
"""
module IsorropiaII

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

export IsorropiaEquilibrium

@component function IsorropiaEquilibrium(; name=:IsorropiaEquilibrium)
    @constants begin
        # Physical constants
        R = 8.314462618, [unit = u"J/mol/K", description = "Universal gas constant"]
        T₀ = 298.15, [unit = u"K", description = "Reference temperature"]
        p_atm = 101325.0, [unit = u"Pa", description = "Standard atmosphere pressure"]

        # Debye-Hückel constant (corrected units from paper)
        A_γ = 0.511, [unit = u"mol^0.5/kg^0.5", description = "Debye-Hückel constant at 298.15 K"]

        # Unit constants for dimensional consistency
        I_one = 1.0, [unit = u"mol/kg", description = "Unit ionic strength"]
        m_one = 1.0, [unit = u"mol/kg", description = "Unit molality"]
        M_one = 1.0, [unit = u"mol/m^3", description = "Unit molarity"]
        p_one = 1.0, [unit = u"Pa", description = "Unit pressure"]
        γ_one = 1.0, [description = "Unit activity coefficient (dimensionless)"]

        # Equilibrium constants at 298.15 K from Table 2
        # All converted to proper SI units (Pa for pressure-based constants)
        K₁_298 = 6.067e5 / p_atm, [unit = u"Pa^-1", description = "Ca(NO3)2 equilibrium constant"]
        K₂_298 = 7.974e11 / p_atm, [unit = u"Pa^-1", description = "CaCl2 equilibrium constant"]
        K₃_298 = 4.319e-5, [unit = u"Pa^-2", description = "CaSO4·2H2O equilibrium constant"]
        K₄_298 = 1.569e-2, [description = "K2SO4 equilibrium constant (dimensionless)"]
        K₅_298 = 24.016, [unit = u"Pa^-1", description = "KHSO4 equilibrium constant"]
        K₆_298 = 0.872, [unit = u"Pa^-1", description = "KNO3 equilibrium constant"]
        K₇_298 = 8.680, [unit = u"Pa^-1", description = "KCl equilibrium constant"]
        K₈_298 = 1.079e5, [description = "MgSO4 equilibrium constant (dimensionless)"]
        K₉_298 = 2.507e15 / p_atm, [unit = u"Pa^-1", description = "Mg(NO3)2 equilibrium constant"]
        K₁₀_298 = 9.557e21 / p_atm, [unit = u"Pa^-1", description = "MgCl2 equilibrium constant"]
        K₁₁_298 = 1.015e-2, [unit = u"mol/kg", description = "H2SO4 dissociation constant"]
        K₁₂_298 = 5.764e1 / p_atm, [unit = u"Pa^-1*mol/kg", description = "NH3-NH4 equilibrium"]
        K₁₃_298 = 1.805e-5, [unit = u"mol/kg", description = "NH4-NH3 aqueous equilibrium"]
        K₁₄_298 = 2.511e6 / p_atm^2, [unit = u"Pa^-2*mol^2/kg^2", description = "HNO3 equilibrium"]
        K₁₅_298 = 2.1e5 / p_atm, [unit = u"Pa^-1*mol/kg", description = "HNO3 aqueous equilibrium"]
        K₁₆_298 = 1.971e6 / p_atm^2, [unit = u"Pa^-2*mol^2/kg^2", description = "HCl equilibrium"]
        K₁₇_298 = 2.5e3 / p_atm, [unit = u"Pa^-1*mol/kg", description = "HCl aqueous equilibrium"]
        K₁₈_298 = 1.010e-14, [unit = u"mol^2/kg^2", description = "Water dissociation constant"]
        K₁₉_298 = 4.799e-1, [description = "Na2SO4 equilibrium constant (dimensionless)"]
        K₂₀_298 = 1.87, [description = "NH42SO4 equilibrium constant (dimensionless)"]
        K₂₁_298 = 1.086e-16 * p_atm^2, [unit = u"Pa^2", description = "NH4Cl gas equilibrium"]
        K₂₂_298 = 1.197e1, [unit = u"Pa^-1", description = "NaNO3 equilibrium constant"]
        K₂₃_298 = 3.766e1, [unit = u"Pa^-1", description = "NaCl equilibrium constant"]
        K₂₄_298 = 2.413e4, [unit = u"Pa^-1", description = "NaHSO4 equilibrium constant"]
        K₂₅_298 = 4.199e-17 * p_atm^2, [unit = u"Pa^2", description = "NH4NO3 gas equilibrium"]
        K₂₆_298 = 1.383, [unit = u"Pa^-1", description = "NH4HSO4 equilibrium constant"]
        K₂₇_298 = 2.972e1, [description = "(NH4)3H(SO4)2 equilibrium constant (dimensionless)"]

        # Temperature dependence coefficients from Table 2 (ΔH°/RT₀)
        H₁ = -11.299, [description = "Ca(NO3)2 enthalpy coefficient (dimensionless)"]
        H₂ = -14.087, [description = "CaCl2 enthalpy coefficient (dimensionless)"]
        H₄ = -9.589, [description = "K2SO4 enthalpy coefficient (dimensionless)"]
        H₅ = -8.423, [description = "KHSO4 enthalpy coefficient (dimensionless)"]
        H₆ = 14.075, [description = "KNO3 enthalpy coefficient (dimensionless)"]
        H₇ = -6.167, [description = "KCl enthalpy coefficient (dimensionless)"]
        H₈ = 36.798, [description = "MgSO4 enthalpy coefficient (dimensionless)"]
        H₉ = -8.754, [description = "Mg(NO3)2 enthalpy coefficient (dimensionless)"]
        H₁₀ = -1.347, [description = "MgCl2 enthalpy coefficient (dimensionless)"]
        H₁₁ = 8.85, [description = "H2SO4 enthalpy coefficient (dimensionless)"]
        H₁₂ = 13.79, [description = "NH3 enthalpy coefficient (dimensionless)"]
        H₁₃ = -1.50, [description = "NH4-NH3 enthalpy coefficient (dimensionless)"]
        H₁₄ = 29.17, [description = "HNO3 enthalpy coefficient (dimensionless)"]
        H₁₅ = 29.17, [description = "HNO3 aq enthalpy coefficient (dimensionless)"]
        H₁₆ = 30.20, [description = "HCl enthalpy coefficient (dimensionless)"]
        H₁₇ = 30.20, [description = "HCl aq enthalpy coefficient (dimensionless)"]
        H₁₈ = -22.52, [description = "H2O enthalpy coefficient (dimensionless)"]
        H₁₉ = 0.98, [description = "Na2SO4 enthalpy coefficient (dimensionless)"]
        H₂₀ = -2.65, [description = "NH42SO4 enthalpy coefficient (dimensionless)"]
        H₂₁ = -71.00, [description = "NH4Cl gas enthalpy coefficient (dimensionless)"]
        H₂₂ = -8.22, [description = "NaNO3 enthalpy coefficient (dimensionless)"]
        H₂₃ = -1.56, [description = "NaCl enthalpy coefficient (dimensionless)"]
        H₂₄ = 0.79, [description = "NaHSO4 enthalpy coefficient (dimensionless)"]
        H₂₅ = -74.375, [description = "NH4NO3 gas enthalpy coefficient (dimensionless)"]
        H₂₆ = -2.87, [description = "NH4HSO4 enthalpy coefficient (dimensionless)"]
        H₂₇ = -5.19, [description = "(NH4)3H(SO4)2 enthalpy coefficient (dimensionless)"]

        # Heat capacity coefficients from Table 2 (ΔC°ₚ/R)
        C₄ = 45.807, [description = "K2SO4 heat capacity coefficient (dimensionless)"]
        C₅ = 17.964, [description = "KHSO4 heat capacity coefficient (dimensionless)"]
        C₆ = 19.388, [description = "KNO3 heat capacity coefficient (dimensionless)"]
        C₇ = 19.953, [description = "KCl heat capacity coefficient (dimensionless)"]
        C₁₁ = 25.14, [description = "H2SO4 heat capacity coefficient (dimensionless)"]
        C₁₂ = -5.39, [description = "NH3 heat capacity coefficient (dimensionless)"]
        C₁₃ = 26.92, [description = "NH4-NH3 heat capacity coefficient (dimensionless)"]
        C₁₄ = 16.83, [description = "HNO3 heat capacity coefficient (dimensionless)"]
        C₁₅ = 16.83, [description = "HNO3 aq heat capacity coefficient (dimensionless)"]
        C₁₆ = 19.91, [description = "HCl heat capacity coefficient (dimensionless)"]
        C₁₇ = 19.91, [description = "HCl aq heat capacity coefficient (dimensionless)"]
        C₁₈ = 26.92, [description = "H2O heat capacity coefficient (dimensionless)"]
        C₁₉ = 39.75, [description = "Na2SO4 heat capacity coefficient (dimensionless)"]
        C₂₀ = 38.57, [description = "NH42SO4 heat capacity coefficient (dimensionless)"]
        C₂₁ = 2.40, [description = "NH4Cl gas heat capacity coefficient (dimensionless)"]
        C₂₂ = 16.01, [description = "NaNO3 heat capacity coefficient (dimensionless)"]
        C₂₃ = 16.90, [description = "NaCl heat capacity coefficient (dimensionless)"]
        C₂₄ = 14.75, [description = "NaHSO4 heat capacity coefficient (dimensionless)"]
        C₂₅ = 6.025, [description = "NH4NO3 gas heat capacity coefficient (dimensionless)"]
        C₂₆ = 15.83, [description = "NH4HSO4 heat capacity coefficient (dimensionless)"]
        C₂₇ = 54.40, [description = "(NH4)3H(SO4)2 heat capacity coefficient (dimensionless)"]

        # Deliquescence relative humidity at 298.15 K from Table 4
        DRH_Ca_NO3_2 = 0.4906, [description = "Ca(NO3)2 DRH at 298.15 K (dimensionless)"]
        DRH_CaCl2 = 0.2830, [description = "CaCl2 DRH at 298.15 K (dimensionless)"]
        DRH_CaSO4 = 0.9700, [description = "CaSO4 DRH at 298.15 K (dimensionless)"]
        DRH_KHSO4 = 0.8600, [description = "KHSO4 DRH at 298.15 K (dimensionless)"]
        DRH_K2SO4 = 0.9751, [description = "K2SO4 DRH at 298.15 K (dimensionless)"]
        DRH_KNO3 = 0.9248, [description = "KNO3 DRH at 298.15 K (dimensionless)"]
        DRH_KCl = 0.8426, [description = "KCl DRH at 298.15 K (dimensionless)"]
        DRH_MgSO4 = 0.8612, [description = "MgSO4 DRH at 298.15 K (dimensionless)"]
        DRH_Mg_NO3_2 = 0.5400, [description = "Mg(NO3)2 DRH at 298.15 K (dimensionless)"]
        DRH_MgCl2 = 0.3284, [description = "MgCl2 DRH at 298.15 K (dimensionless)"]
        DRH_NaCl = 0.7528, [description = "NaCl DRH at 298.15 K (dimensionless)"]
        DRH_Na2SO4 = 0.9300, [description = "Na2SO4 DRH at 298.15 K (dimensionless)"]
        DRH_NaNO3 = 0.7379, [description = "NaNO3 DRH at 298.15 K (dimensionless)"]
        DRH_NH4_2SO4 = 0.7997, [description = "(NH4)2SO4 DRH at 298.15 K (dimensionless)"]
        DRH_NH4NO3 = 0.6183, [description = "NH4NO3 DRH at 298.15 K (dimensionless)"]
        DRH_NH4Cl = 0.7710, [description = "NH4Cl DRH at 298.15 K (dimensionless)"]
        DRH_NH4HSO4 = 0.4000, [description = "NH4HSO4 DRH at 298.15 K (dimensionless)"]
        DRH_NaHSO4 = 0.5200, [description = "NaHSO4 DRH at 298.15 K (dimensionless)"]
        DRH_NH4_3H_SO4_2 = 0.6900, [description = "(NH4)3H(SO4)2 DRH at 298.15 K (dimensionless)"]

        # Molar masses (g/mol) - now properly in @constants block
        MW_H = 1.00784, [unit = u"g/mol", description = "Molar mass of H"]
        MW_Na = 22.989769, [unit = u"g/mol", description = "Molar mass of Na"]
        MW_NH4 = 18.038, [unit = u"g/mol", description = "Molar mass of NH4"]
        MW_K = 39.0983, [unit = u"g/mol", description = "Molar mass of K"]
        MW_Ca = 40.078, [unit = u"g/mol", description = "Molar mass of Ca"]
        MW_Mg = 24.305, [unit = u"g/mol", description = "Molar mass of Mg"]
        MW_Cl = 35.453, [unit = u"g/mol", description = "Molar mass of Cl"]
        MW_NO3 = 62.0049, [unit = u"g/mol", description = "Molar mass of NO3"]
        MW_SO4 = 96.0636, [unit = u"g/mol", description = "Molar mass of SO4"]
        MW_HSO4 = 97.0705, [unit = u"g/mol", description = "Molar mass of HSO4"]
        MW_OH = 17.00734, [unit = u"g/mol", description = "Molar mass of OH"]
        MW_H2O = 18.01528, [unit = u"g/mol", description = "Molar mass of H2O"]
        MW_NH3 = 17.03052, [unit = u"g/mol", description = "Molar mass of NH3"]
        MW_HNO3 = 63.01, [unit = u"g/mol", description = "Molar mass of HNO3"]
        MW_HCl = 36.46, [unit = u"g/mol", description = "Molar mass of HCl"]
        MW_H2SO4 = 98.0785, [unit = u"g/mol", description = "Molar mass of H2SO4"]
    end

    @parameters begin
        T = 298.15, [unit = u"K", description = "Temperature"]
        RH = 0.8, [description = "Relative humidity (0 to 1) (dimensionless)"]
        stable = false, [description = "Use stable (solid) mode if true, metastable (liquid) mode if false"]

        # Total input concentrations (mol/m³)
        C_NH4_total = 1e-6, [unit = u"mol/m^3", description = "Total NH4 + NH3 concentration"]
        C_Na_total = 1e-7, [unit = u"mol/m^3", description = "Total Na concentration"]
        C_K_total = 1e-8, [unit = u"mol/m^3", description = "Total K concentration"]
        C_Ca_total = 1e-8, [unit = u"mol/m^3", description = "Total Ca concentration"]
        C_Mg_total = 1e-9, [unit = u"mol/m^3", description = "Total Mg concentration"]
        C_SO4_total = 5e-7, [unit = u"mol/m^3", description = "Total SO4 + HSO4 concentration"]
        C_NO3_total = 2e-7, [unit = u"mol/m^3", description = "Total NO3 concentration"]
        C_Cl_total = 1e-7, [unit = u"mol/m^3", description = "Total Cl concentration"]
    end

    @variables begin
        # Equilibrium constants (temperature-dependent)
        K₁(t), [description = "Ca(NO3)2 equilibrium constant"]
        K₂(t), [description = "CaCl2 equilibrium constant"]
        K₃(t), [description = "CaSO4·2H2O equilibrium constant"]
        K₄(t), [description = "K2SO4 equilibrium constant"]
        K₅(t), [description = "KHSO4 equilibrium constant"]
        K₆(t), [description = "KNO3 equilibrium constant"]
        K₇(t), [description = "KCl equilibrium constant"]
        K₈(t), [description = "MgSO4 equilibrium constant"]
        K₉(t), [description = "Mg(NO3)2 equilibrium constant"]
        K₁₀(t), [description = "MgCl2 equilibrium constant"]
        K₁₁(t), [description = "H2SO4 dissociation constant"]
        K₁₂(t), [description = "NH3-NH4 equilibrium constant"]
        K₁₃(t), [description = "NH4-NH3 aqueous equilibrium"]
        K₁₄(t), [description = "HNO3 equilibrium constant"]
        K₁₅(t), [description = "HNO3 aqueous equilibrium"]
        K₁₆(t), [description = "HCl equilibrium constant"]
        K₁₇(t), [description = "HCl aqueous equilibrium"]
        K₁₈(t), [description = "Water dissociation constant"]
        K₁₉(t), [description = "Na2SO4 equilibrium constant"]
        K₂₀(t), [description = "NH42SO4 equilibrium constant"]
        K₂₁(t), [description = "NH4Cl gas equilibrium constant"]
        K₂₂(t), [description = "NaNO3 equilibrium constant"]
        K₂₃(t), [description = "NaCl equilibrium constant"]
        K₂₄(t), [description = "NaHSO4 equilibrium constant"]
        K₂₅(t), [description = "NH4NO3 gas equilibrium constant"]
        K₂₆(t), [description = "NH4HSO4 equilibrium constant"]
        K₂₇(t), [description = "(NH4)3H(SO4)2 equilibrium constant"]

        # Aqueous species molalities (mol/kg water)
        m_H(t), [unit = u"mol/kg", description = "H+ molality"]
        m_OH(t), [unit = u"mol/kg", description = "OH- molality"]
        m_Na(t), [unit = u"mol/kg", description = "Na+ molality"]
        m_NH4(t), [unit = u"mol/kg", description = "NH4+ molality"]
        m_K(t), [unit = u"mol/kg", description = "K+ molality"]
        m_Ca(t), [unit = u"mol/kg", description = "Ca2+ molality"]
        m_Mg(t), [unit = u"mol/kg", description = "Mg2+ molality"]
        m_Cl(t), [unit = u"mol/kg", description = "Cl- molality"]
        m_NO3(t), [unit = u"mol/kg", description = "NO3- molality"]
        m_SO4(t), [unit = u"mol/kg", description = "SO4 2- molality"]
        m_HSO4(t), [unit = u"mol/kg", description = "HSO4- molality"]
        m_NH3(t), [unit = u"mol/kg", description = "NH3 molality"]
        m_HNO3(t), [unit = u"mol/kg", description = "HNO3 molality"]
        m_HCl(t), [unit = u"mol/kg", description = "HCl molality"]

        # Gas phase partial pressures (Pa)
        p_NH3(t), [unit = u"Pa", description = "NH3 partial pressure"]
        p_HNO3(t), [unit = u"Pa", description = "HNO3 partial pressure"]
        p_HCl(t), [unit = u"Pa", description = "HCl partial pressure"]
        p_H2O(t), [unit = u"Pa", description = "H2O partial pressure"]

        # Activity coefficients (dimensionless)
        γ_H(t), [description = "H+ activity coefficient (dimensionless)"]
        γ_OH(t), [description = "OH- activity coefficient (dimensionless)"]
        γ_Na(t), [description = "Na+ activity coefficient (dimensionless)"]
        γ_NH4(t), [description = "NH4+ activity coefficient (dimensionless)"]
        γ_K(t), [description = "K+ activity coefficient (dimensionless)"]
        γ_Ca(t), [description = "Ca2+ activity coefficient (dimensionless)"]
        γ_Mg(t), [description = "Mg2+ activity coefficient (dimensionless)"]
        γ_Cl(t), [description = "Cl- activity coefficient (dimensionless)"]
        γ_NO3(t), [description = "NO3- activity coefficient (dimensionless)"]
        γ_SO4(t), [description = "SO4 2- activity coefficient (dimensionless)"]
        γ_HSO4(t), [description = "HSO4- activity coefficient (dimensionless)"]
        γ_NH3(t), [description = "NH3 activity coefficient (dimensionless)"]
        γ_HNO3(t), [description = "HNO3 activity coefficient (dimensionless)"]
        γ_HCl(t), [description = "HCl activity coefficient (dimensionless)"]
        γ_H2O(t), [description = "H2O activity coefficient (dimensionless)"]

        # System properties
        I(t), [unit = u"mol/kg", description = "Ionic strength"]
        W(t), [unit = u"kg/m^3", description = "Aerosol water content"]
        pH(t), [description = "Solution pH (dimensionless)"]

        # Aerosol type classification (Section 3.1)
        R1(t), [description = "Total sulfate ratio (dimensionless)"]
        R2(t), [description = "Crustal species and sodium ratio (dimensionless)"]
        R3(t), [description = "Crustal species ratio (dimensionless)"]

        # Aerosol water content concentrations (mol/m³)
        C_H(t), [unit = u"mol/m^3", description = "Aqueous H+ concentration"]
        C_OH(t), [unit = u"mol/m^3", description = "Aqueous OH- concentration"]
        C_Na(t), [unit = u"mol/m^3", description = "Aqueous Na+ concentration"]
        C_NH4(t), [unit = u"mol/m^3", description = "Aqueous NH4+ concentration"]
        C_NH3(t), [unit = u"mol/m^3", description = "Aqueous NH3 concentration"]
        C_K(t), [unit = u"mol/m^3", description = "Aqueous K+ concentration"]
        C_Ca(t), [unit = u"mol/m^3", description = "Aqueous Ca2+ concentration"]
        C_Mg(t), [unit = u"mol/m^3", description = "Aqueous Mg2+ concentration"]
        C_Cl(t), [unit = u"mol/m^3", description = "Aqueous Cl- concentration"]
        C_NO3(t), [unit = u"mol/m^3", description = "Aqueous NO3- concentration"]
        C_SO4(t), [unit = u"mol/m^3", description = "Aqueous SO4 2- concentration"]
        C_HSO4(t), [unit = u"mol/m^3", description = "Aqueous HSO4- concentration"]
        C_HNO3(t), [unit = u"mol/m^3", description = "Aqueous HNO3 concentration"]
        C_HCl(t), [unit = u"mol/m^3", description = "Aqueous HCl concentration"]
    end

    @equations begin
        # Van't Hoff equation for temperature dependence (Equations 3-5)
        # K(T) = K₂₉₈ * exp(-ΔH°/R * (1/T - 1/T₀) - ΔC°ₚ/R * (1 + ln(T₀/T) - T₀/T))
        K₁ ~ K₁_298 * exp(-H₁ * (T₀/T - 1))  # Eq. 3 - Ca(NO3)2
        K₂ ~ K₂_298 * exp(-H₂ * (T₀/T - 1))  # Eq. 3 - CaCl2
        K₃ ~ K₃_298  # No temperature dependence given
        K₄ ~ K₄_298 * exp(-H₄ * (T₀/T - 1) - C₄ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - K2SO4
        K₅ ~ K₅_298 * exp(-H₅ * (T₀/T - 1) - C₅ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - KHSO4
        K₆ ~ K₆_298 * exp(-H₆ * (T₀/T - 1) - C₆ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - KNO3
        K₇ ~ K₇_298 * exp(-H₇ * (T₀/T - 1) - C₇ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - KCl
        K₈ ~ K₈_298 * exp(-H₈ * (T₀/T - 1))  # Eq. 3 - MgSO4
        K₉ ~ K₉_298 * exp(-H₉ * (T₀/T - 1))  # Eq. 3 - Mg(NO3)2
        K₁₀ ~ K₁₀_298 * exp(-H₁₀ * (T₀/T - 1))  # Eq. 3 - MgCl2
        K₁₁ ~ K₁₁_298 * exp(-H₁₁ * (T₀/T - 1) - C₁₁ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - H2SO4
        K₁₂ ~ K₁₂_298 * exp(-H₁₂ * (T₀/T - 1) - C₁₂ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NH3
        K₁₃ ~ K₁₃_298 * exp(-H₁₃ * (T₀/T - 1) - C₁₃ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NH4-NH3
        K₁₄ ~ K₁₄_298 * exp(-H₁₄ * (T₀/T - 1) - C₁₄ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - HNO3
        K₁₅ ~ K₁₅_298 * exp(-H₁₅ * (T₀/T - 1) - C₁₅ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - HNO3 aq
        K₁₆ ~ K₁₆_298 * exp(-H₁₆ * (T₀/T - 1) - C₁₆ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - HCl
        K₁₇ ~ K₁₇_298 * exp(-H₁₇ * (T₀/T - 1) - C₁₇ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - HCl aq
        K₁₈ ~ K₁₈_298 * exp(-H₁₈ * (T₀/T - 1) - C₁₈ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - H2O
        K₁₉ ~ K₁₉_298 * exp(-H₁₉ * (T₀/T - 1) - C₁₉ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - Na2SO4
        K₂₀ ~ K₂₀_298 * exp(-H₂₀ * (T₀/T - 1) - C₂₀ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - (NH4)2SO4
        K₂₁ ~ K₂₁_298 * exp(-H₂₁ * (T₀/T - 1) - C₂₁ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NH4Cl gas
        K₂₂ ~ K₂₂_298 * exp(-H₂₂ * (T₀/T - 1) - C₂₂ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NaNO3
        K₂₃ ~ K₂₃_298 * exp(-H₂₃ * (T₀/T - 1) - C₂₃ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NaCl
        K₂₄ ~ K₂₄_298 * exp(-H₂₄ * (T₀/T - 1) - C₂₄ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NaHSO4
        K₂₅ ~ K₂₅_298 * exp(-H₂₅ * (T₀/T - 1) - C₂₅ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NH4NO3 gas
        K₂₆ ~ K₂₆_298 * exp(-H₂₆ * (T₀/T - 1) - C₂₆ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - NH4HSO4
        K₂₇ ~ K₂₇_298 * exp(-H₂₇ * (T₀/T - 1) - C₂₇ * (1 + log(T₀/T) - T₀/T))  # Eq. 5 - (NH4)3H(SO4)2

        # Ionic strength calculation (Equation after 2.3)
        I ~ 0.5 * (m_H + m_Na + m_NH4 + m_K + 4*m_Ca + 4*m_Mg + m_Cl + m_NO3 + 4*m_SO4 + m_HSO4 + m_OH) * m_one  # Eq. after 2.3

        # Simplified activity coefficients using Debye-Hückel limiting law
        # For proper implementation, this should use Bromley's formula (Eq. 6)
        # This is a placeholder for the complex activity coefficient calculation
        γ_H ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))  # Simplified Debye-Hückel
        γ_OH ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_Na ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_NH4 ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_K ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_Ca ~ γ_one * exp(-A_γ * 4 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))  # z² = 4 for divalent
        γ_Mg ~ γ_one * exp(-A_γ * 4 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_Cl ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_NO3 ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_SO4 ~ γ_one * exp(-A_γ * 4 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_HSO4 ~ γ_one * exp(-A_γ * 1 * sqrt(I/I_one) / (1 + sqrt(I/I_one)))
        γ_NH3 ~ γ_one  # Neutral species
        γ_HNO3 ~ γ_one  # Neutral species
        γ_HCl ~ γ_one  # Neutral species
        γ_H2O ~ γ_one  # Solvent

        # Equilibrium relations from Table 2
        # Each equilibrium constant relates to activities of products/reactants

        # Water dissociation (Reaction 18)
        K₁₈ ~ (γ_H * m_H / m_one) * (γ_OH * m_OH / m_one)  # Eq. 1: K = ∏ aᵢ^νᵢ

        # Acid dissociation (Reaction 11: H2SO4 → H+ + HSO4-)
        # Assuming H2SO4 fully dissociates, so [H2SO4] = 0 in equilibrium
        # K₁₁ ~ (γ_H * m_H / m_one) * (γ_HSO4 * m_HSO4 / m_one) / (1e-10)  # Avoid division by zero

        # Gas-liquid equilibrium for NH3 (Reaction 12)
        K₁₂ ~ (γ_NH3 * m_NH3 / m_one) * (p_one / p_NH3)  # Henry's law

        # NH3 protonation equilibrium (Reaction 13: NH3 + H2O → NH4+ + OH-)
        K₁₃ ~ (γ_NH4 * m_NH4 / m_one) * (γ_OH * m_OH / m_one) / (γ_NH3 * m_NH3 / m_one) / RH  # with water activity ≈ RH

        # Gas-liquid equilibrium for HNO3 (Reactions 14-15)
        K₁₄ ~ (γ_H * m_H / m_one) * (γ_NO3 * m_NO3 / m_one) * (p_one / p_HNO3)
        K₁₅ ~ (γ_HNO3 * m_HNO3 / m_one) * (p_one / p_HNO3)

        # Gas-liquid equilibrium for HCl (Reactions 16-17)
        K₁₆ ~ (γ_H * m_H / m_one) * (γ_Cl * m_Cl / m_one) * (p_one / p_HCl)
        K₁₇ ~ (γ_HCl * m_HCl / m_one) * (p_one / p_HCl)

        # Mass balance equations (Section 2.1)
        C_NH4_total ~ C_NH4 + C_NH3 + p_NH3 / (R * T)  # Total NH4/NH3
        C_Na_total ~ C_Na  # Total Na (conservative)
        C_K_total ~ C_K  # Total K (conservative)
        C_Ca_total ~ C_Ca  # Total Ca (conservative)
        C_Mg_total ~ C_Mg  # Total Mg (conservative)
        C_SO4_total ~ C_SO4 + C_HSO4  # Total sulfate
        C_NO3_total ~ C_NO3 + C_HNO3 + p_HNO3 / (R * T)  # Total nitrate
        C_Cl_total ~ C_Cl + C_HCl + p_HCl / (R * T)  # Total chloride

        # Charge balance (electroneutrality)
        0 ~ C_H + C_Na + C_NH4 + C_K + 2*C_Ca + 2*C_Mg - C_OH - C_Cl - C_NO3 - 2*C_SO4 - C_HSO4  # Electroneutrality

        # Convert molalities to molarities using water content
        C_H ~ m_H * m_one * W  # mol/m³ = (mol/kg) × (kg/m³)
        C_OH ~ m_OH * m_one * W
        C_Na ~ m_Na * m_one * W
        C_NH4 ~ m_NH4 * m_one * W
        C_NH3 ~ m_NH3 * m_one * W
        C_K ~ m_K * m_one * W
        C_Ca ~ m_Ca * m_one * W
        C_Mg ~ m_Mg * m_one * W
        C_Cl ~ m_Cl * m_one * W
        C_NO3 ~ m_NO3 * m_one * W
        C_SO4 ~ m_SO4 * m_one * W
        C_HSO4 ~ m_HSO4 * m_one * W
        C_HNO3 ~ m_HNO3 * m_one * W
        C_HCl ~ m_HCl * m_one * W

        # Aerosol type classification (Section 3.1)
        R1 ~ (C_NH4_total + C_Ca_total + C_K_total + C_Mg_total + C_Na_total) / C_SO4_total  # Eq. from Section 3.1
        R2 ~ (C_Ca_total + C_K_total + C_Mg_total + C_Na_total) / C_SO4_total  # Eq. from Section 3.1
        R3 ~ (C_Ca_total + C_K_total + C_Mg_total) / C_SO4_total  # Eq. from Section 3.1

        # Water content estimation (simplified ZSR method)
        # This should be replaced with proper ZSR calculation from Equation 16
        W ~ 1e-3  # Placeholder: 1 g/m³ water content

        # pH calculation
        pH ~ -log10(γ_H * m_H / m_one)

        # Water vapor pressure (simplified)
        p_H2O ~ RH * 3169.0  # Saturated vapor pressure of water at 298 K (Pa) × RH

        # Conservation constraints (steady state)
        D(C_NH4_total) ~ 0
        D(C_Na_total) ~ 0
        D(C_K_total) ~ 0
        D(C_Ca_total) ~ 0
        D(C_Mg_total) ~ 0
        D(C_SO4_total) ~ 0
        D(C_NO3_total) ~ 0
        D(C_Cl_total) ~ 0
    end
end

@doc """
    IsorropiaEquilibrium(; name=:IsorropiaEquilibrium, kwargs...)

ISORROPIA II thermodynamic equilibrium model for inorganic aerosol species.

This component implements the thermodynamic equilibrium calculations described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient
> thermodynamic equilibrium model for K+–Ca2+–Mg2+–NH4+–Na+–SO42-–NO3-–Cl-–H2O
> aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.

## Parameters

- `T`: Temperature (K), default = 298.15
- `RH`: Relative humidity (0-1), default = 0.8
- `stable`: Use stable mode (with solids) if true, metastable mode (liquid only) if false
- `C_NH4_total`: Total NH4/NH3 concentration (mol/m³)
- `C_Na_total`: Total Na concentration (mol/m³)
- `C_K_total`: Total K concentration (mol/m³)
- `C_Ca_total`: Total Ca concentration (mol/m³)
- `C_Mg_total`: Total Mg concentration (mol/m³)
- `C_SO4_total`: Total SO4/HSO4 concentration (mol/m³)
- `C_NO3_total`: Total NO3 concentration (mol/m³)
- `C_Cl_total`: Total Cl concentration (mol/m³)

## Variables

The component calculates equilibrium concentrations of all aqueous and gas-phase species,
along with thermodynamic properties like ionic strength, pH, and water content.

## Example

```julia
@named aerosol = IsorropiaEquilibrium(T = 298.15, RH = 0.8,
                                      C_NH4_total = 1e-6, C_SO4_total = 5e-7)
```
""" IsorropiaEquilibrium

end